#ifndef BATCHER_X_PARTICLE_EMITTER_LOGIC_HPP
#define BATCHER_X_PARTICLE_EMITTER_LOGIC_HPP

#include <unordered_map>
#include "../../graphics/draw_info/draw_info.hpp"
#include "../../graphics/batcher/generated/batcher.hpp"
#include "../../graphics/texture_packer/texture_packer.hpp"
#include "../../utility/print_utils/print_utils.hpp"
#include "../../graphics/particle_emitter/particle_emitter.hpp"
#include "../../graphics/fps_camera/fps_camera.hpp"

class BatcherXParticleEmitterLogic {
  public:
    TexturePacker &texture_packer;
    BoundedUniqueIDGenerator &ltw_object_id_generator;
    Batcher &batcher;

    draw_info::TransformedIVPTPGroup smoke_tig_for_copying;

    BatcherXParticleEmitterLogic(TexturePacker &texture_packer, BoundedUniqueIDGenerator &ltw_object_id_generator,
                                 Batcher &batcher)
        : texture_packer(texture_packer), ltw_object_id_generator(ltw_object_id_generator), batcher(batcher) {

        // TODO: need to use gfp here?
        auto smoke_filepath = "assets/smoke_64px.png";
        auto smoke_st = texture_packer.get_packed_texture_sub_texture(smoke_filepath);

        auto vertices = generate_square_vertices(0, 0, 0.5);
        auto indices = generate_rectangle_indices();

        draw_info::IVPTexturePacked smoke_ivptp(
            indices, vertices, smoke_st.texture_coordinates, smoke_st.texture_coordinates,
            smoke_st.packed_texture_index, smoke_st.packed_texture_bounding_box_index, smoke_filepath,
            batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.object_id_generator.get_id());

        smoke_tig_for_copying = draw_info::TransformedIVPTPGroup({smoke_ivptp}, -1);
    }

    std::unordered_map<int, draw_info::TransformedIVPTPGroup> particle_id_to_tig;
    void on_particle_death(int particle_id) {

        p(std::format("deleting particle: {}", particle_id));
        auto particle_tig = particle_id_to_tig[particle_id];
        ltw_object_id_generator.reclaim_id(particle_tig.id);

        for (const auto &ivptp : particle_tig.ivptps) {
            batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.delete_object(ivptp.id);
        }

        particle_id_to_tig.erase(particle_id);
    }

    void on_particle_spawn(int particle_id) {
        p(std::format("spawning particle: {}", particle_id));
        auto smoke_tig_copy = smoke_tig_for_copying;
        smoke_tig_copy.regenerate_ids(
            ltw_object_id_generator,
            batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.object_id_generator);

        particle_id_to_tig[particle_id] = smoke_tig_copy;
    }

    void draw_particles(ParticleEmitter &particle_emitter, FPSCamera &fps_camera, double delta_time,
                        glm::mat4 projection, glm::mat4 view) {
        particle_emitter.update(delta_time, projection * view);
        auto particles = particle_emitter.get_particles_sorted_by_distance();

        for (size_t i = 0; i < particles.size(); ++i) {
            auto &curr_particle = particles[i];

            //  compute the up vector (assuming we want it to be along the y-axis)
            glm::vec3 up = glm::vec3(0.0f, 1.0f, 0.0f);
            glm::vec3 forward = fps_camera.transform.compute_forward_vector();

            glm::vec3 right = glm::normalize(glm::cross(up, forward));

            up = glm::normalize(glm::cross(forward, right));

            // TODO: replace with functions in transfomr
            // this makes it billboarded
            glm::mat4 rotation_matrix = glm::mat4(1.0f);
            rotation_matrix[0] = glm::vec4(right, 0.0f);
            rotation_matrix[1] = glm::vec4(up, 0.0f);
            rotation_matrix[2] = glm::vec4(-forward, 0.0f); // We negate the direction for correct facing

            glm::mat4 transform = glm::translate(glm::mat4(1.0f), curr_particle.transform.position);
            transform *= rotation_matrix;
            transform = glm::scale(transform, curr_particle.transform.scale);

            auto tig = particle_id_to_tig[curr_particle.id];

            batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.ltw_matrices[tig.id] = transform;

            for (const auto &ivptp : tig.ivptps) {
                std::vector<unsigned int> ltw_mat_idxs(ivptp.xyz_positions.size(), tig.id);
                std::vector<int> ptis(ivptp.xyz_positions.size(), ivptp.packed_texture_index);
                std::vector<int> ptbbis(ivptp.xyz_positions.size(), ivptp.packed_texture_bounding_box_index);
                batcher.texture_packer_cwl_v_transformation_ubos_1024_shader_batcher.queue_draw(
                    ivptp.id, ivptp.indices, ltw_mat_idxs, ptis, ivptp.packed_texture_coordinates, ptbbis,
                    ivptp.xyz_positions, false);
            }
        }
    }
};

#endif // BATCHER_X_PARTICLE_EMITTER_LOGIC_HPP
