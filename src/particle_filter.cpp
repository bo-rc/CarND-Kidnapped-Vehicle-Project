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

#include "particle_filter.h"

using namespace std;

void ParticleFilter::init(double x, double y, double theta, double std[], int num) {
    // TODO: Set the number of particles. Initialize all particles to first position (based on estimates of
    //   x, y, theta and their uncertainties from GPS) and all weights to 1.
    // Add random Gaussian noise to each particle.
    // NOTE: Consult particle_filter.h for more information about this method (and others in this file).
    num_particles = num;

    std::random_device rd{};
    std::mt19937 gen{rd()};

    normal_distribution<double> x_gen(x, std[0]);
    normal_distribution<double> y_gen(y, std[1]);
    normal_distribution<double> theta_gen(theta, std[2]);

    for (int i = 0; i < num_particles; ++i) {
	Particle p;
	p.id = i;
	p.x = x_gen(gen);
	p.y = y_gen(gen);
	p.theta = theta_gen(gen);
	p.weight = 1.0;
	particles.push_back(p);
    }

    is_initialized = true;
}


void ParticleFilter::prediction(double delta_t, double std_pos[], double velocity, double yaw_rate) {

    for (auto& particle: particles) {
	#define EPSILON 0.00001
	if (fabs(yaw_rate) < EPSILON) {
	    particle.x += velocity * delta_t * cos(particle.theta);
	    particle.y += velocity * delta_t * sin(particle.theta);
	} else {
	    particle.x += velocity / yaw_rate * (sin(particle.theta + yaw_rate*delta_t) - sin(particle.theta));
	    particle.y += velocity / yaw_rate * (cos(particle.theta) - cos(particle.theta + yaw_rate*delta_t));
	    particle.theta += yaw_rate * delta_t;
	}

	std::random_device rd{};
	std::mt19937 gen{rd()};

	normal_distribution<double> x_gen(0, std_pos[0]);
	normal_distribution<double> y_gen(0, std_pos[1]);
	normal_distribution<double> theta_gen(0, std_pos[2]);
	particle.x += x_gen(gen);
	particle.y += y_gen(gen);
	particle.theta += theta_gen(gen);
    }
}

void ParticleFilter::dataAssociation(const std::vector<LandmarkObs>& predicted, std::vector<LandmarkObs>& observations) {
    for (auto& ob: observations) {
	int id = -1;
	auto min_distance_squared = numeric_limits<double>::max();

	for (auto& pred: predicted) {
	    auto distance_squared = pow(ob.x - pred.x, 2) + pow(ob.y - pred.y, 2);

	    // relexation algorithm
	    if (distance_squared < min_distance_squared) {
		min_distance_squared = distance_squared;
		id = pred.id;
	    }
	}
	ob.id = id;
    }
}

// sensor_range = 50m, std_landmark: x:0.3, y:0.3, observations are noisy, map_landmarks are known
void ParticleFilter::updateWeights(double sensor_range, double std_landmark[], 
				   const std::vector<LandmarkObs> &observations, const Map &map_landmarks) {
    /*
     * For each particle
     */
    for (auto& particle: particles) {
	/*
	 * finding all landmarks from map that are in range
	 * this is the predicted landmark set for this particle
	 */
	vector<LandmarkObs> landmarks_in_range_of_particle;
	for (const auto& lm: map_landmarks.landmark_list) {
	    if (pow(lm.x_f - particle.x, 2) + pow(lm.y_f - particle.y, 2) < pow(sensor_range,2))
		landmarks_in_range_of_particle.emplace_back(lm.id_i, lm.x_f, lm.y_f);
	}

	/*
	 * tranform observations to map coordinates
	 */
	vector<LandmarkObs> landmarks_observed;
	for (const auto& ob: observations) {
	    auto x = cos(particle.theta)*ob.x - sin(particle.theta)*ob.y + particle.x;
	    auto y = sin(particle.theta)*ob.x + cos(particle.theta)*ob.y + particle.y;
	    landmarks_observed.emplace_back(ob.id, x, y);
	}

	/*
	 * data association
	 */
	dataAssociation(landmarks_in_range_of_particle, landmarks_observed);

	/*
	 * calculate weight for this particle
	 */
	particle.weight = 1.0;
	for (const auto& landmark_observed: landmarks_observed) {
	    auto predicted_x = 0.0;
	    auto predicted_y = 0.0;

	    // get corresponding predicted landmark of particle
	    for (const auto& landmark_particle: landmarks_in_range_of_particle) {
		if (landmark_particle.id == landmark_observed.id) {
		    predicted_x = landmark_particle.x;
		    predicted_y = landmark_particle.y;
		    break;
		}
	    }

	    // weight for this observation
	    // https://en.wikipedia.org/wiki/Gaussian_function
	    auto weight = (1/(2*M_PI*std_landmark[0]*std_landmark[1])) * exp(-(pow(predicted_x-landmark_observed.x,2)/(2*pow(std_landmark[0],2))+(pow(predicted_y-landmark_observed.y,2)/(2*pow(std_landmark[1],2)))));
	    particle.weight *= weight;
	}
    }
}

void ParticleFilter::resample() {
    /*
     * weights list
     */
    vector<double> weights;
    for (int i = 0; i < num_particles; ++i) {
	weights.push_back(particles[i].weight);
    }

    random_device rd;
    mt19937 gen(rd());
    discrete_distribution<> d(begin(weights), end(weights));

    /*
     * draw the same number of particles
     */
    vector<Particle> new_particles;
    for (int i = 0; i < num_particles; ++i) {
	auto idx = d(gen);
	new_particles.push_back(particles[idx]);
    }

    swap(particles, new_particles);
}

Particle ParticleFilter::SetAssociations(Particle& particle, const std::vector<int>& associations, 
					 const std::vector<double>& sense_x, const std::vector<double>& sense_y)
{
    //particle: the particle to assign each listed association, and association's (x,y) world coordinates mapping to
    // associations: The landmark id that goes along with each listed association
    // sense_x: the associations x mapping already converted to world coordinates
    // sense_y: the associations y mapping already converted to world coordinates

    particle.associations= associations;
    particle.sense_x = sense_x;
    particle.sense_y = sense_y;
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
