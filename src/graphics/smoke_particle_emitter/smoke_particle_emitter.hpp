#ifndef SMOKE_PARTICLE_EMITTER_HPP
#define SMOKE_PARTICLE_EMITTER_HPP

#include <functional>
#include <glm/glm.hpp>
#include <random>
#include "../particle_emitter/particle_emitter.hpp"

class SmokeParticleEmitter {
  public:
    ParticleEmitter particle_emitter;

    SmokeParticleEmitter(std::function<void(int)> on_particle_spawn_callback,
                         std::function<void(int)> on_particle_death_callback)
        : particle_emitter(life_span_lambda(), initial_velocity_lambda(), velocity_change_lambda(), scaling_lambda(),
                           rotation_lambda(), spawn_delay_lambda(), on_particle_spawn_callback,
                           on_particle_death_callback) {}

  private:
    static std::function<float()> life_span_lambda() {
        return []() -> float {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> dist(1.0f, 3.0f);
            return dist(rng);
        };
    }

    static std::function<glm::vec3()> initial_velocity_lambda() {
        return []() -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.5f, 0.5f); // Small lateral variance
            std::uniform_real_distribution<float> upward_dist(1.0f, 2.0f);      // Strong upward push

            // Initial upward push with slight lateral drift
            float dx = horizontal_dist(rng);
            float dy = upward_dist(rng); // Main upward velocity
            float dz = horizontal_dist(rng);

            return glm::vec3(dx, dy, dz);
            /*return glm::vec3(0, 0, 0);*/
        };
    }

    static std::function<glm::vec3(float, float)> velocity_change_lambda() {
        return [](float life_percentage, float delta_time) -> glm::vec3 {
            static std::mt19937 rng(std::random_device{}());
            std::uniform_real_distribution<float> horizontal_dist(-0.010f, 0.010f); // Small lateral variance
            std::uniform_real_distribution<float> vertical_dist(0.2f, 0.3f);

            float accel_x = horizontal_dist(rng);
            float accel_y = vertical_dist(rng);
            float accel_z = horizontal_dist(rng);

            glm::vec3 smoke_push_down = -glm::vec3(accel_x, accel_y, accel_z) * delta_time;
            return smoke_push_down;
        };
    }

    static std::function<float(float)> scaling_lambda() {
        return [](float life_percentage) -> float { return std::max(life_percentage * 2.0f, 0.0f); };
    }

    static std::function<float(float)> rotation_lambda() {
        return [](float life_percentage) -> float { return life_percentage / 5.0f; };
    }

    static std::function<float()> spawn_delay_lambda() {
        return []() -> float {
            return 0.01f; // Spawn a new particle every 0.05 seconds
        };
    }
};

#endif // SMOKE_PARTICLE_EMITTER_HPP
