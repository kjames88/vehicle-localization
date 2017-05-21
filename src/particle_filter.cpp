/*
 * particle_filter.cpp
 *
 *  Created on: Dec 12, 2016
 *      Author: Tiffany Huang
 */

#include <random>
#include <algorithm>
#include <iostream>
#include <numeric>
#include <math.h> 
#include <iostream>
#include <sstream>
#include <string>
#include <iterator>
#include <assert.h>

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
	// TODO: Set the number of particles. Initialize all particles to first position (based on estimates of 
	//   x, y, theta and their uncertainties from GPS) and all weights to 1. 
	// Add random Gaussian noise to each particle.
	// NOTE: Consult particle_filter.h for more information about this method (and others in this file).
  default_random_engine gen;
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);
  num_particles = 1000;
  for (int i=0; i<num_particles; i++) {
    Particle p;
    p.id = i;
    p.x = dist_x(gen);
    p.y = dist_y(gen);
    p.theta = dist_theta(gen);
    p.weight = 1.0;
    
    particles.push_back(p);
  }
  
  is_initialized = true;


  // std::vector<std::vector<int> > test;
  // std::vector<int> l1 = {1, 2};
  // std::vector<int> l2 = {3, 6, 5};
  // std::vector<int> l3 = {6, 7};
  // std::vector<int> l4 = {8, 9, 10, 11};
  // test.push_back(l1);
  // test.push_back(l2);
  // test.push_back(l3);
  // test.push_back(l4);
  // test = recursive_assoc(test, 0);
  // for (auto it=test.begin(); it!=test.end(); it++) {
  //   for (auto it2=it->begin(); it2!=it->end(); it2++)
  //     std::cout << *it2 << ", ";
  //   std::cout << endl;
  // }
}

void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {
	// TODO: Add measurements to each particle and add random Gaussian noise.
	// NOTE: When adding noise you may find std::normal_distribution and std::default_random_engine useful.
	//  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
	//  http://www.cplusplus.com/reference/random/default_random_engine/
  default_random_engine gen;
  normal_distribution<double> dist_v(velocity, std_pos[0]);    // FIXME
  normal_distribution<double> dist_yaw(yaw_rate, std_pos[1]);  // FIXME
  double vtd = velocity / yaw_rate;
  double tdt = delta_t * yaw_rate;
  for (int i=0; i<num_particles; i++) {
    double theta_0 = particles.at(i).theta;
    double x_f = particles.at(i).x + (vtd * (sin(theta_0 + tdt) - sin(theta_0)));
    double y_f = particles.at(i).y + (vtd * (cos(theta_0) - cos(theta_0 + tdt)));
    double theta_f = theta_0 + tdt;
    particles.at(i).x = x_f;
    particles.at(i).y = y_f;
    particles.at(i).theta = theta_f;
  }
}

void ParticleFilter::dataAssociation(std::vector<LandmarkObs> predicted, std::vector<LandmarkObs>& observations) {
	// TODO: Find the predicted measurement that is closest to each observed measurement and assign the 
	//   observed measurement to this particular landmark.
	// NOTE: this method will NOT be called by the grading code. But you will probably find it useful to 
	//   implement this method and use it as a helper during the updateWeights phase.

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
		std::vector<LandmarkObs> observations, Map map_landmarks) {
	// TODO: Update the weights of each particle using a mult-variate Gaussian distribution. You can read
	//   more about this distribution here: https://en.wikipedia.org/wiki/Multivariate_normal_distribution
	// NOTE: The observations are given in the VEHICLE'S coordinate system. Your particles are located
	//   according to the MAP'S coordinate system. You will need to transform between the two systems.
	//   Keep in mind that this transformation requires both rotation AND translation (but no scaling).
	//   The following is a good resource for the theory:
	//   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
	//   and the following is a good resource for the actual equation to implement (look at equation 
	//   3.33. Note that you'll need to switch the minus sign in that equation to a plus to account 
	//   for the fact that the map's y-axis actually points downwards.)
	//   http://planning.cs.uiuc.edu/node99.html
  
  double sigma[2][2] = {std_landmark[0], 0, 0, std_landmark[1]};
  int m = observations.size();
  
  for (int p=0; p<num_particles; p++) {
    Particle particle = particles.at(p);
    double xt = particle.x;
    double yt = particle.y;
    std::vector<LandmarkObs> obs_world;
    for (auto it = observations.begin(); it != observations.end(); it++) {
      // particle has world coordinates but observation is relative to particle coordinates
      // translate observation to world coordinates
      //   rotation followed by translation
      // Rotate by (desired - actual) = (0 - particle theta)
      // Counter-clockwise rotation with inverted y:  [cos(a) sin(a), -sin(a) cos(a)]
      double xobs = it->x;
      double yobs = it->y;
      double theta = 0 - particle.theta;  // 0 degrees is facing upward (14.13 key)
      double xw = (xobs * cos(theta) + yobs * sin(theta)) + xt;
      double yw = (-1.0 * xobs * sin(theta) + yobs * cos(theta)) + yt;
      LandmarkObs obs;
      obs.x = xw;
      obs.y = yw;
      obs.id = 0;
      obs_world.push_back(obs);

      std::cout << "particle x=" << xt << " y=" << yt << " theta=" << particle.theta <<
        " obs x=" << xobs << " y=" << yobs << " world obs x=" << xw << " y=" << yw << std::endl;
    }
    
    // Establish correspondences
    //   Filter reasonable pointwise correspondence to limit the number of permutations
    //   Find the maximum weight among permutations and assign to the particle
    std::vector<std::vector<int> > candidates;
    candidates.resize(m);
    int c = 0;
    for (auto ito = obs_world.begin(); ito != obs_world.end(); ito++) {
      // Exclude any observations that are out of sensor range
      double dist = sqrt(pow(ito->x - xt, 2) + pow(ito->y - yt, 2));      
      if (dist <= sensor_range) {
        for (auto itm = map_landmarks.landmark_list.begin(); itm != map_landmarks.landmark_list.end(); itm++) {
          // Exclude any landmarks that are out of sensor range
          double dist = sqrt(pow(itm->x_f - xt, 2) + pow(itm->y_f - yt, 2));
          int landmark_id = itm->id_i;
          if (dist <= sensor_range) {
            // Compute the landmark and observation angles and exclude any landmark if angle differs by +/- pi/8 radians
            double rho = atan2(ito->y, ito->x);
            double landmark_rho = atan2(itm->y_f, itm->x_f);
            if (abs(rho - landmark_rho) < M_PI/8) {
              candidates.at(c).push_back(landmark_id);  // list of landmark IDs that could associate with this observation
              std::cout << "accept landmark " << landmark_id << " world obs x=" << ito->x << " y=" << ito->y << std::endl;
            }
          }
        }
      }
      c++;
    }
    // Compute multivariate gaussian for each combination of reasonable associations
    // Take the maximum product as the correspondence and weight for the particle
    // In case of a missing landmark set the weight to a minimum value

    bool fail = false;
    for (int i=0; i<m; i++) {
      if (candidates.at(i).empty())
        fail = true;
    }
    if (fail == false) {
      double w_max = 0.0;
      std::vector<std::vector<int> > sets = recursive_assoc(candidates, 0);
      for (auto it=sets.begin(); it!=sets.end(); it++) {
        double w = mv_gaussian(obs_world, *it, sigma, map_landmarks);
        if (w > w_max) {
          w_max = w;
          particles.at(p).associations = *it;
        }
      }
      weights[p] = w_max;
    } else {
      weights[p] = 0.0001;
    }
    particles.at(p).weight = weights[p];
  }
  
}

std::vector<std::vector<int> > ParticleFilter::recursive_assoc(std::vector<std::vector<int> > candidates, int idx) {
  std::vector<vector<int> > r;
  if (idx == candidates.size() - 1) {
    for (auto it=candidates.at(idx).begin(); it != candidates.at(idx).end(); it++) {
      std::vector<int> e;
      e.push_back(*it);
      r.push_back(e);
    }
  } else {
    std::vector<vector<int> > below = recursive_assoc(candidates, idx + 1);
    for (int i=0; i<candidates.at(idx).size(); i++) {
      for (auto it = below.begin(); it != below.end(); it++) {
        // don't create any sets with duplicate entries
        bool proceed = true;
        for (auto it2 = it->begin(); it2 != it->end(); it2++) {
          if (*it2 == candidates.at(idx).at(i))
            proceed = false;
        }
        if (proceed) {
          std::vector<int> extended;
          extended.push_back(candidates.at(idx).at(i));
          extended.insert(extended.end(), it->begin(), it->end());
          r.push_back(extended);
        }
      }
    }
  }
  return r;
}

double ParticleFilter::mv_gaussian(std::vector<LandmarkObs> const& observations,
                                   std::vector<int> const& associations,
                                   double sigma[2][2],
                                   Map const& landmark_map) {

  // 2x2 matrix inverse:
  //   http://www.mathcentre.ac.uk/resources/uploaded/sigma-matrices7-2009-1.pdf

  double mult = 1.0 / (sigma[0][0] * sigma[1][1] - sigma[0][1] * sigma[1][0]);
  double sigma_inv[2][2] = {mult * sigma[1][1], -1.0 * mult * sigma[0][1], -1.0 * mult * sigma[1][0], mult * sigma[0][0]};
  
  assert(observations.size() == associations.size());
  double w = 0.0;
  int m = observations.size();
  for (int i=0; i<m; i++) {
    double delta[2];
    Map::single_landmark_s landmark = landmark_map.landmark_list.at(associations.at(i));
    delta[0] = observations.at(i).x - landmark.x_f;
    delta[1] = observations.at(i).y - landmark.y_f;
    double e1[2] = {sigma_inv[0][0] * delta[0] + sigma_inv[0][1] * delta[1],
                    sigma_inv[1][0] * delta[0] + sigma_inv[1][1] * delta[1]};
    double e = -0.5 * (delta[0] * e1[0] + delta[1] * e1[1]);
    e = exp(e);
  }
  
  return w;
}

void ParticleFilter::resample() {
	// TODO: Resample particles with replacement with probability proportional to their weight. 
	// NOTE: You may find std::discrete_distribution helpful here.
	//   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution

}

Particle ParticleFilter::SetAssociations(Particle particle, std::vector<int> associations, std::vector<double> sense_x, std::vector<double> sense_y)
{
	//particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
	// associations: The landmark id that goes along with each listed association
	// sense_x: the associations x mapping already converted to world coordinates
	// sense_y: the associations y mapping already converted to world coordinates

	//Clear the previous associations
	particle.associations.clear();
	particle.sense_x.clear();
	particle.sense_y.clear();

	particle.associations= associations;
 	particle.sense_x = sense_x;
 	particle.sense_y = sense_y;

 	return particle;
}

string ParticleFilter::getAssociations(Particle best)
{
	vector<int> v = best.associations;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<int>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseX(Particle best)
{
	vector<double> v = best.sense_x;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
string ParticleFilter::getSenseY(Particle best)
{
	vector<double> v = best.sense_y;
	stringstream ss;
    copy( v.begin(), v.end(), ostream_iterator<float>(ss, " "));
    string s = ss.str();
    s = s.substr(0, s.length()-1);  // get rid of the trailing space
    return s;
}
