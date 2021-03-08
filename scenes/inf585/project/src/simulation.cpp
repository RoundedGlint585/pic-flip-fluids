//
// Created by roundedglint585 on 3/3/21.
//

#include "simulation.h"

void simulate(Particles particles, float dt)
{

	// Update values
	particles.step(dt);

	// Numerical integration
	/*float const damping = 0.005f;
	size_t const N = particles.size();
	float const m = sph_parameters.m;
	for (size_t k = 0; k < N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;
		vec3& f = particles[k].f;

		v = (1 - damping) * v + dt * f / m;
		p = p + dt * v;
	}


	// Collision
	float const epsilon = 1e-3f;
	for (size_t k = 0; k < N; ++k)
	{
		vec3& p = particles[k].p;
		vec3& v = particles[k].v;

		// small perturbation to avoid alignment
		if (p.y < -1) { p.y = -1 + epsilon * rand_interval();  v.y *= -0.5f; }
		if (p.x < -1) { p.x = -1 + epsilon * rand_interval();  v.x *= -0.5f; }
		if (p.x > 1) { p.x = 1 - epsilon * rand_interval();  v.x *= -0.5f; }
	}*/

}