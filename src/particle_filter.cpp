/**
 * particle_filter.cpp
 *
 * Created on: Dec 12, 2016
 * Author: Tiffany Huang
 */

#include "particle_filter.h"

#include <math.h>
#include <algorithm>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <vector>

#include "helper_functions.h"

using std::string;
using std::vector;
using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[]) {
  /**
   * TODO: Set the number of particles. Initialize all particles to 
   *   first position (based on estimates of x, y, theta and their uncertainties
   *   from GPS) and all weights to 1. 
   * TODO: Add random Gaussian noise to each particle.
   * NOTE: Consult particle_filter.h for more information about this method 
   *   (and others in this file).
   */

  default_random_engine gen;

  //create normal distribution for x, y and theta
  normal_distribution<double> dist_x(x, std[0]);
  normal_distribution<double> dist_y(y, std[1]);
  normal_distribution<double> dist_theta(theta, std[2]);


  num_particles = 100;  // TODO: Set the number of particles

  // Sample and form these normal distrubtions like this: 
  // x = dist_x(gen);
  // where "gen" is the random engine initialized earlier.

  for (int i=0; i<num_particles; ++i) {
      Particle p;
      p.id = i;
      p.x = dist_x(gen);
      p.y = dist_y(gen);
      p.theta = dist_theta(gen);
      p.weight = 1.0;
      particles.push_back(p);
    }

    is_initialized = true; //forgot to write this the first time

}

void ParticleFilter::prediction(double delta_t, double std_pos[], 
                                double velocity, double yaw_rate) {
  /**
   * TODO: Add measurements to each particle and add random Gaussian noise.
   * NOTE: When adding noise you may find std::normal_distribution 
   *   and std::default_random_engine useful.
   *  http://en.cppreference.com/w/cpp/numeric/random/normal_distribution
   *  http://www.cplusplus.com/reference/random/default_random_engine/
   */

  default_random_engine gen;

  normal_distribution<double> noise_x(0, std_pos[0]);
  normal_distribution<double> noise_y(0, std_pos[1]);
  normal_distribution<double> noise_theta(0, std_pos[2]);


  for (int i=0; i<num_particles; ++i) {
    if( fabs(yaw_rate) < 0.0001){  // constant velocity, did not specify this situation the first time
      particles[i].x = particles[i].x + velocity * delta_t * cos(particles[i].theta);
      particles[i].y = particles[i].y + velocity * delta_t * sin(particles[i].theta);

    } else {
      particles[i].x = particles[i].x + velocity / yaw_rate * (sin(particles[i].theta + yaw_rate * delta_t) - sin(particles[i].theta)) + noise_x(gen);
      particles[i].y = particles[i].y + velocity / yaw_rate * (cos(particles[i].theta) - cos(particles[i].theta + yaw_rate * delta_t));
      particles[i].theta = particles[i].theta + yaw_rate * delta_t;
    }

    //add sensor noise
    particles[i].x = particles[i].x + noise_x(gen);
    particles[i].y = particles[i].y + noise_y(gen);
    particles[i].theta = particles[i].theta + noise_theta(gen);
  }

}

void ParticleFilter::dataAssociation(vector<LandmarkObs> predicted, 
                                     vector<LandmarkObs>& observations) {
  /**
   * TODO: Find the predicted measurement that is closest to each 
   *   observed measurement and assign the observed measurement to this 
   *   particular landmark.
   * NOTE: this method will NOT be called by the grading code. But you will 
   *   probably find it useful to implement this method and use it as a helper 
   *   during the updateWeights phase.
   */
  for (int i=0; i<observations.size(); ++i) {
       double min_dist = numeric_limits<double>::max();
    for (int j=0; j<predicted.size(); ++j) {
       double distance = dist(predicted[j].x, predicted[j].y, observations[i].x, observations[i].y);
        if (distance < min_dist) {
          min_dist = distance;
          observations[i].id = predicted[j].id;
        }
    }
  }

}

void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
                                   const vector<LandmarkObs> &observations, 
                                   const Map &map_landmarks) {
  /**
   * TODO: Update the weights of each particle using a mult-variate Gaussian 
   *   distribution. You can read more about this distribution here: 
   *   https://en.wikipedia.org/wiki/Multivariate_normal_distribution
   * NOTE: The observations are given in the VEHICLE'S coordinate system. 
   *   Your particles are located according to the MAP'S coordinate system. 
   *   You will need to transform between the two systems. Keep in mind that
   *   this transformation requires both rotation AND translation (but no scaling).
   *   The following is a good resource for the theory:
   *   https://www.willamette.edu/~gorr/classes/GeneralGraphics/Transforms/transforms2d.htm
   *   and the following is a good resource for the actual equation to implement
   *   (look at equation 3.33) http://planning.cs.uiuc.edu/node99.html
   */

  for (int k=0; k<num_particles; ++k) {

    particles[k].weight = 1.0;

    vector<LandmarkObs> predictions; //landmarks within sensor range of the particle
    vector<LandmarkObs> observations_os; //list of observations transformed from vehicle cooridinates to map coordinates

    //pick landmarks just within sensor range
    for (int i=0; i<map_landmarks.landmark_list.size(); ++i) {
      int lm_id = map_landmarks.landmark_list[i].id_i;
      double lm_x = map_landmarks.landmark_list[i].x_f;
      double lm_y = map_landmarks.landmark_list[i].y_f;
      double distance = dist(lm_x, lm_y, particles[k].x, particles[k].y);
      if (distance < sensor_range) {
        predictions.push_back(LandmarkObs{lm_id, lm_x, lm_y});
      }
    }

  //transform observations in the car coordinates system to map coordinates
    for (int i=0; i<observations.size(); ++i) {
      double x_t = particles[k].x + (cos(particles[k].theta) * observations[i].x) - (sin(particles[k].theta) * observations[i].y);
      double y_t = particles[k].y + (sin(particles[k].theta) * observations[i].x) + (cos(particles[k].theta) * observations[i].y);
      observations_os.push_back(LandmarkObs{observations[i].id, x_t, y_t});
    }

  // find landmark index for each observation
    dataAssociation(predictions, observations_os);

  // compute particle's weight
    double prediction_x;
    double prediction_y;

    for (int i=0; i<observations_os.size(); ++i) {
      for (int j=0; j<predictions.size(); ++j) {
        if (predictions[j].id == observations_os[i].id) {
          prediction_x = predictions[j].x;
          prediction_y = predictions[j].y;
        }
      }
    double weight_obs = 1.0 / (2.0 * M_PI * std_landmark[0] * std_landmark[1]) * exp (-((observations_os[i].x - prediction_x) * (observations_os[i].x - prediction_x) / (2.0 * std_landmark[0] * std_landmark [0]) + \
                    (observations_os[i].y - prediction_y) * (observations_os[i].y - prediction_y) / (2.0 * std_landmark[1] * std_landmark [1])));
    particles[k].weight = particles[k].weight * weight_obs;
    }
    weights.push_back(particles[k].weight);
  }

}

void ParticleFilter::resample() {
  /**
   * TODO: Resample particles with replacement with probability proportional 
   *   to their weight. 
   * NOTE: You may find std::discrete_distribution helpful here.
   *   http://en.cppreference.com/w/cpp/numeric/random/discrete_distribution
   */
    default_random_engine gen;
    discrete_distribution<> dist(begin(weights), end(weights));

    // create resampled particles
    vector<Particle> resampled_particles;
    resampled_particles.resize(num_particles);

    // resample the particles according to weights
    for(int i=0; i<num_particles; i++){
      int idx = dist(gen);
      resampled_particles[i] = particles[idx];
    }

    // assign the resampled_particles to the previous particles
    particles = resampled_particles;

    // clear the weight vector for the next round
    weights.clear();
}

void ParticleFilter::SetAssociations(Particle& particle, 
                                     const vector<int>& associations, 
                                     const vector<double>& sense_x, 
                                     const vector<double>& sense_y) {
  // particle: the particle to which assign each listed association, 
  //   and association's (x,y) world coordinates mapping
  // associations: The landmark id that goes along with each listed association
  // sense_x: the associations x mapping already converted to world coordinates
  // sense_y: the associations y mapping already converted to world coordinates
  particle.associations= associations;
  particle.sense_x = sense_x;
  particle.sense_y = sense_y;
}

string ParticleFilter::getAssociations(Particle best) {
  vector<int> v = best.associations;
  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<int>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}

string ParticleFilter::getSenseCoord(Particle best, string coord) {
  vector<double> v;

  if (coord == "X") {
    v = best.sense_x;
  } else {
    v = best.sense_y;
  }

  std::stringstream ss;
  copy(v.begin(), v.end(), std::ostream_iterator<float>(ss, " "));
  string s = ss.str();
  s = s.substr(0, s.length()-1);  // get rid of the trailing space
  return s;
}