/*
 * Copyright (C) 2019, unclearness
 * All rights reserved.
 */

#include "currender/rasterizer.h"

#include <array>
#include <cassert>
#include <unordered_set>

#include "currender/pixel_shader.h"
#include "currender/timer.h"
#include "currender/util_private.h"

namespace {

template <typename T>
void Argsort(const std::vector<T>& data, std::vector<size_t>* indices) {
  indices->resize(data.size());
  std::iota(indices->begin(), indices->end(), 0);
  std::sort(indices->begin(), indices->end(),
            [&data](size_t i1, size_t i2) { return data[i1] < data[i2]; });
}

inline float EdgeFunction(const Eigen::Vector3f& a, const Eigen::Vector3f& b,
                          const Eigen::Vector3f& c) {
  return (c[0] - a[0]) * (b[1] - a[1]) - (c[1] - a[1]) * (b[0] - a[0]);
}

inline bool InsideTri(float x, float y){
  Eigen::Vector3f pixel_sample(x, y,
                               0.0f);
  float w0 = EdgeFunction(v1_i, v2_i, pixel_sample);
  float w1 = EdgeFunction(v2_i, v0_i, pixel_sample);
  float w2 = EdgeFunction(v0_i, v1_i, pixel_sample);


 if ((w0 >= 0 && w1 >= 0 && w2 >= 0) ||
      (w0 <= 0 && w1 <= 0 && w2 <= 0)) {
    return true;
  }

  return false;

}


inline void InitMinMaxTableNaive(
    const Eigen::Vector3f& v0_i, const Eigen::Vector3f& v1_i,
    const Eigen::Vector3f& v2_i, const Eigen::Vector3f& face_normal,
    uint32_t x0, uint32_t x1, uint32_t y0, uint32_t y1,
    std::shared_ptr<const currender::Camera> camera,
    std::vector<std::pair<int, int>>* minmax_table) {
  for (int y = y0; y <= y1; ++y) {
    std::pair<int, int>& table_row = (*minmax_table)[y - y0];
    table_row.first = x1+1;
    table_row.second = x0-1;

    for (int x = x0; x <= x1; ++x) {
      Eigen::Vector3f pixel_sample(static_cast<float>(x), static_cast<float>(y),
                                   0.0f);
      float w0 = EdgeFunction(v1_i, v2_i, pixel_sample);
      float w1 = EdgeFunction(v2_i, v0_i, pixel_sample);
      float w2 = EdgeFunction(v0_i, v1_i, pixel_sample);

      Eigen::Vector3f ray_w;
      camera->ray_w(static_cast<int>(x), static_cast<int>(y), &ray_w);
      // even if back-face culling is enabled, dont' skip back-face
      // need to update z-buffer to handle front-face occluded by back-face
      bool backface = face_normal.dot(ray_w) > 0;

      if ((!backface && (w0 >= 0 && w1 >= 0 && w2 >= 0)) ||
          (backface && (w0 <= 0 && w1 <= 0 && w2 <= 0))) {
        table_row.first = std::min(table_row.first, x);
        table_row.second = std::max(table_row.second, x);
      }
    }
  }
}

inline void InitMinMaxTable(const Eigen::Vector3f& v0_i,
                            const Eigen::Vector3f& v1_i,
                            const Eigen::Vector3f& v2_i,
                            const Eigen::Vector3f& face_normal, uint32_t x0,
                            uint32_t x1, uint32_t y0, uint32_t y1,
                            std::shared_ptr<const currender::Camera> camera,
                            std::vector<std::pair<int, int>>* minmax_table) {
  std::array<Eigen::Vector3f, 3> vertices{{v0_i, v1_i, v2_i}};
  std::unordered_set<int> indices{{0, 1, 2}};
  int ymin_index = static_cast<int>(std::distance(
      vertices.begin(),
      std::min_element(vertices.begin(), vertices.end(),
                       [](const Eigen::Vector3f& a, const Eigen::Vector3f& b) {
                         return a.y() < b.y();
                       })));
  const Eigen::Vector3f& v_ymin = vertices[ymin_index];
  indices.erase(ymin_index);

  int l_index = ymin_index != 0 ? 0 : 1;
  int r_index = l_index;
  Eigen::Vector3f v_l{vertices[l_index]}, v_r{vertices[r_index]};
  for (const int& i : indices) {
    const Eigen::Vector3f& v = vertices[i];
    if (v.x() < v_l.x()) {
      v_l = v;
      l_index = i;
    }
    if (v_r.x() < v.x()) {
      v_r = v;
      r_index = i;
    }
  }

  //     v_y_min *
  //            /  _
  //           /     _
  //          /        _
  //         *
  //

#if 0
				  currender::LOGI("v_ymin (%f, %f, %f)\n", v_ymin.x(), v_ymin.y(), v_ymin.z());
  currender::LOGI("v_l (%f, %f, %f)\n", v_l.x(), v_l.y(), v_l.z());
  currender::LOGI("v_r (%f, %f, %f)\n", v_r.x(), v_r.y(), v_r.z());

#endif  // 0

  for (int y = y0; y <= y1; ++y) {
    std::pair<int, int>& table_row = (*minmax_table)[y - y0];
    table_row.first = x1 + 1;
    table_row.second = x0 - 1;
  }

  // v_ymin to v_l line
  // update x_min
  {
    Eigen::Vector3f delta = v_l - v_ymin;
    double grad = double(delta.x()) / double(delta.y());
    //assert(grad > 0);
    for (int y = y0; y <= static_cast<int>(std::floor(v_l.y())); ++y) {
      std::pair<int, int>& table_row = (*minmax_table)[y - y0];
      double xf = double(grad) * double(y- v_ymin.y()) + double(v_ymin.x());
#if 0
				      float integral_part{0.0f};
      float fractional_part = std::modf(xf, &integral_part);
      int xi = static_cast<int>(integral_part);
      if (fractional_part < 0.5f) {
        table_row.first = xi + 1;
      } else {
        table_row.first = xi + 1;
      }
#else
      int xi = static_cast<int>(std::ceil(xf));
      table_row.first = xi;
#endif  // 0

    }
  }

  // v_y_min to v_r line
  // update x_max
  {
    Eigen::Vector3f delta = v_r - v_ymin;
    double grad = double(delta.x()) / double(delta.y());
   // assert(grad > 0);
    for (int y = y0; y <= static_cast<int>(std::floor(v_r.y())); ++y) {
      std::pair<int, int>& table_row = (*minmax_table)[y - y0];
      double xf = double(grad) * double(y - v_ymin.y()) + double(v_ymin.x());
#if 0
      float integral_part{0.0f};
      float fractional_part = std::modf(xf, &integral_part);
      int xi = static_cast<int>(integral_part);
								      if (fractional_part < 0.5f) {
        table_row.second = xi;
      } else {
        table_row.second = xi - 1;
      }
#else
      int xi = static_cast<int>(std::floor(xf));
      table_row.second = xi;
#endif  // 1

      //table_row.second = xi - 1;
    }
  }

  // final line: v_l to v_r
  if (v_l.y() < v_r.y()) {
    // v_l is upper than v_r case
    // update x_min
    Eigen::Vector3f delta = v_r - v_l;
    float grad = delta.x() / delta.y();

    int y_start =
        std::max(static_cast<int>(y0), static_cast<int>(std::ceil(v_l.y())));
    int y_end = std::min(static_cast<int>(y1), static_cast<int>(std::floor(v_r.y())));

    //assert(grad > 0);
    for (int y = y_start;
         y <= y_end; ++y) {
      std::pair<int, int>& table_row = (*minmax_table)[y - y0];
      float xf = grad * (y - v_l.y()) + v_l.x();
#if 0
				      float integral_part{0.0f};
      float fractional_part = std::modf(xf, &integral_part);
      int xi = static_cast<int>(integral_part);
      if (fractional_part < 0.5f) {
        table_row.first = xi + 1;
      } else {
        table_row.first = xi + 1;
      }
#else
      int xi = static_cast<int>(std::ceil(xf));
      table_row.first = xi;
#endif  // 0

    }

  } else {
    // v_l is lower than v_r case
    // update x_max
    Eigen::Vector3f delta = v_l - v_r;
    float grad = delta.x() / delta.y();
   // assert(grad < 0);
    int y_start =
        std::max(static_cast<int>(y0), static_cast<int>(std::ceil(v_r.y())));
    int y_end =
        std::min(static_cast<int>(y1), static_cast<int>(std::floor(v_l.y())));

    for (int y = y_start;
         y <= y_end; ++y) {
      std::pair<int, int>& table_row = (*minmax_table)[y - y0];
      float xf = grad * (y -v_r.y())+ v_r.x();
#if 0
      float integral_part{0.0f};
      float fractional_part = std::modf(xf, &integral_part);
      int xi = static_cast<int>(integral_part);
								      if (fractional_part < 0.5f) {
        table_row.second = xi;
      } else {
        table_row.second = xi - 1;
      }
#else
      int xi = static_cast<int>(std::floor(xf));
      table_row.second = xi;
#endif  // 

    }
  }

  for (int y = y0; y <= y1; ++y) {
    std::pair<int, int>& table_row = (*minmax_table)[y - y0];
    table_row.first = std::max(static_cast<int>(x0), table_row.first);
    table_row.second = std::min(static_cast<int>(x1), table_row.second);
  }

}

}  // namespace

namespace currender {

// Rasterizer::Impl implementation
class Rasterizer::Impl {
  bool mesh_initialized_{false};
  std::shared_ptr<const Camera> camera_{nullptr};
  std::shared_ptr<const Mesh> mesh_{nullptr};
  RendererOption option_;

 public:
  Impl();
  ~Impl();

  explicit Impl(const RendererOption& option);
  void set_option(const RendererOption& option);

  void set_mesh(std::shared_ptr<const Mesh> mesh);

  bool PrepareMesh();

  void set_camera(std::shared_ptr<const Camera> camera);

  bool Render(Image3b* color, Image1f* depth, Image3f* normal, Image1b* mask,
              Image1i* face_id) const;

  bool RenderColor(Image3b* color) const;
  bool RenderDepth(Image1f* depth) const;
  bool RenderNormal(Image3f* normal) const;
  bool RenderMask(Image1b* mask) const;
  bool RenderFaceId(Image1i* face_id) const;

  bool RenderW(Image3b* color, Image1w* depth, Image3f* normal, Image1b* mask,
               Image1i* face_id) const;
  bool RenderDepthW(Image1w* depth) const;
};

Rasterizer::Impl::Impl() {}
Rasterizer::Impl::~Impl() {}

Rasterizer::Impl::Impl(const RendererOption& option) { set_option(option); }

void Rasterizer::Impl::set_option(const RendererOption& option) {
  option.CopyTo(&option_);
}

void Rasterizer::Impl::set_mesh(std::shared_ptr<const Mesh> mesh) {
  mesh_initialized_ = false;
  mesh_ = mesh;

  if (mesh_->face_normals().empty()) {
    LOGW("face normal is empty. culling and shading may not work\n");
  }

  if (mesh_->normals().empty()) {
    LOGW("vertex normal is empty. shading may not work\n");
  }
}

bool Rasterizer::Impl::PrepareMesh() {
  if (mesh_ == nullptr) {
    LOGE("mesh has not been set\n");
    return false;
  }

  mesh_initialized_ = true;

  return true;
}

void Rasterizer::Impl::set_camera(std::shared_ptr<const Camera> camera) {
  camera_ = camera;
}

bool Rasterizer::Impl::Render(Image3b* color, Image1f* depth, Image3f* normal,
                              Image1b* mask, Image1i* face_id) const {
  if (!ValidateAndInitBeforeRender(mesh_initialized_, camera_, mesh_, option_,
                                   color, depth, normal, mask, face_id)) {
    return false;
  }

  // make pixel shader
  std::unique_ptr<PixelShader> pixel_shader = PixelShaderFactory::Create(
      option_.diffuse_color, option_.interp, option_.diffuse_shading);

  OrenNayarParam oren_nayar_param(option_.oren_nayar_sigma);

  const Eigen::Matrix3f w2c_R = camera_->w2c().rotation().cast<float>();
  const Eigen::Vector3f w2c_t = camera_->w2c().translation().cast<float>();

  Timer<> timer;
  timer.Start();

  // project face to 2d (fully parallel)
  std::vector<Eigen::Vector3f> camera_vertices(mesh_->vertices().size());
  std::vector<Eigen::Vector3f> camera_normals(mesh_->vertices().size());
  std::vector<float> camera_depth_list(mesh_->vertices().size());
  std::vector<Eigen::Vector3f> image_vertices(mesh_->vertices().size());

  // get projected vertex positions
  for (int i = 0; i < static_cast<int>(mesh_->vertices().size()); i++) {
    camera_vertices[i] = w2c_R * mesh_->vertices()[i] + w2c_t;
    camera_normals[i] = w2c_R * mesh_->normals()[i];
    camera_depth_list[i] = camera_vertices[i].z();
    camera_->Project(camera_vertices[i], &image_vertices[i]);
  }

  Image1f depth_internal;
  Image1f* depth_{depth};
  if (depth_ == nullptr) {
    depth_ = &depth_internal;
  }
  depth_->Init(camera_->width(), camera_->height(), 0.0f);

  Image1i face_id_internal;
  Image1i* face_id_{face_id};
  if (face_id_ == nullptr) {
    face_id_ = &face_id_internal;
  }
  face_id_->Init(camera_->width(), camera_->height(), -1);

  // 255: backface, 0:frontface
  Image1b backface_image(camera_->width(), camera_->height(), 0);

  // 0:(1 - u - v), 1:u, 2:v
  Image3f weight_image(camera_->width(), camera_->height(), 0.0f);

  // init once with enough size
  std::vector<std::pair<int, int>> minmax_table(camera_->height());

  // make face id image by z-buffer method
  for (int i = 0; i < static_cast<int>(mesh_->vertex_indices().size()); i++) {
    const Eigen::Vector3i& face = mesh_->vertex_indices()[i];
    const Eigen::Vector3f& v0_i = image_vertices[face[0]];
    const Eigen::Vector3f& v1_i = image_vertices[face[1]];
    const Eigen::Vector3f& v2_i = image_vertices[face[2]];

    // skip if a vertex is back of the camera
    // todo: add near and far plane
    if (v0_i.z() < 0.0f || v1_i.z() < 0.0f || v2_i.z() < 0.0f) {
      continue;
    }

    float xmin = std::min({v0_i.x(), v1_i.x(), v2_i.x()});
    float ymin = std::min({v0_i.y(), v1_i.y(), v2_i.y()});
    float xmax = std::max({v0_i.x(), v1_i.x(), v2_i.x()});
    float ymax = std::max({v0_i.y(), v1_i.y(), v2_i.y()});

    // the triangle is out of screen
    if (xmin > camera_->width() - 1 || xmax < 0 ||
        ymin > camera_->height() - 1 || ymax < 0) {
      continue;
    }

    uint32_t x0 = std::max(int32_t(0), (int32_t)(std::ceil(xmin)));
    uint32_t x1 = std::min(camera_->width() - 1, (int32_t)(std::floor(xmax)));
    uint32_t y0 = std::max(int32_t(0), (int32_t)(std::ceil(ymin)));
    uint32_t y1 = std::min(camera_->height() - 1, (int32_t)(std::floor(ymax)));

    float area = EdgeFunction(v0_i, v1_i, v2_i);
    if (std::abs(area) < std::numeric_limits<float>::min()) {
      continue;
    }
    float inv_denom = 1.0f / area;

    const auto& face_normal = mesh_->face_normals()[i];
    if (i == 27723) {
      LOGI("id %d\n", i);

#if 0

    InitMinMaxTableNaive(v0_i, v1_i, v2_i, face_normal, x0, x1, y0, y1, camera_,
                         &minmax_table);
#endif
#if 0
				    LOGI("naive\n");
    for (uint32_t y = y0; y <= y1; ++y) {
      const std::pair<int, int>& table_row = minmax_table[y - y0];
      LOGI("%d %d (%d, %d)\n", y - y0, y, table_row.first, table_row.second);
    }
#endif
    }
#if 1
    InitMinMaxTable(v0_i, v1_i, v2_i, face_normal, x0, x1, y0, y1, camera_,
                         &minmax_table);
#endif
#if 0
    LOGI("new\n");
    for (uint32_t y = y0; y <= y1; ++y) {
      const std::pair<int, int>& table_row = minmax_table[y - y0];
      LOGI("%d %d (%d, %d)\n", y - y0, y, table_row.first, table_row.second);
    }
    LOGI("\n");

#endif  // 1

    for (uint32_t y = y0; y <= y1; ++y) {
      const std::pair<int, int>& table_row = minmax_table[y - y0];
      for (uint32_t x = table_row.first; x <= table_row.second; ++x) {
        Eigen::Vector3f ray_w;
        camera_->ray_w(static_cast<int>(x), static_cast<int>(y), &ray_w);
        // even if back-face culling is enabled, dont' skip back-face
        // need to update z-buffer to handle front-face occluded by back-face
        bool backface = face_normal.dot(ray_w) > 0;

        Eigen::Vector3f pixel_sample(static_cast<float>(x),
                                     static_cast<float>(y), 0.0f);
        float u = inv_denom * EdgeFunction(v2_i, v0_i, pixel_sample);
        float v = inv_denom * EdgeFunction(v0_i, v1_i, pixel_sample);

        //u = std::min(std::max(0.0f, u), 1.0f);
        //v = std::min(std::max(0.0f, v), 1.0f);
        assert(u >= 0 && u <= 1.0);
        assert(v >= 0 && v <= 1.0);
#if 0
         //original
          pixel_sample.z() = w0 * v0_i.z() + w1 * v1_i.z() + w2 * v2_i.z();
#else
        /** Perspective-Correct Interpolation **/
        float depthnorm_w = (1.0f - u - v) / v0_i.z();
        float depthnorm_u = u / v1_i.z();
        float depthnorm_v = v / v2_i.z();

        pixel_sample.z() = 1.0f / (depthnorm_u + depthnorm_v + depthnorm_w);

        float correct_w = depthnorm_w * pixel_sample.z();
        float correct_u = depthnorm_u * pixel_sample.z();
        float correct_v = depthnorm_v * pixel_sample.z();
        /** Perspective-Correct Interpolation **/
#endif
        if (depth_->at(x, y, 0) < std::numeric_limits<float>::min() ||
            pixel_sample.z() < depth_->at(x, y, 0)) {
          depth_->at(x, y, 0) = pixel_sample.z();
          face_id_->at(x, y, 0) = i;
          weight_image.at(x, y, 0) = correct_w;
          weight_image.at(x, y, 1) = correct_u;
          weight_image.at(x, y, 2) = correct_v;
          backface_image.at(x, y, 0) = backface ? 255 : 0;
        }
      }
    }
  }

  // make images by referring to face id image
  for (int y = 0; y < backface_image.height(); y++) {
    for (int x = 0; x < backface_image.width(); x++) {
      if (option_.backface_culling && backface_image.at(x, y, 0) == 255) {
        depth_->at(x, y, 0) = 0.0f;
        face_id_->at(x, y, 0) = -1;
        continue;
      }

      int fid = face_id_->at(x, y, 0);
      if (fid > 0) {
        Eigen::Vector3f ray_w;
        camera_->ray_w(x, y, &ray_w);

        float w0 = weight_image.at(x, y, 0);
        float w1 = weight_image.at(x, y, 1);
        float w2 = weight_image.at(x, y, 2);

        // fill mask
        if (mask != nullptr) {
          mask->at(x, y, 0) = 255;
        }

        // calculate shading normal
        Eigen::Vector3f shading_normal_w{0.0f, 0.0f, 0.0f};
        if (option_.shading_normal == ShadingNormal::kFace) {
          shading_normal_w = mesh_->face_normals()[fid];
        } else if (option_.shading_normal == ShadingNormal::kVertex) {
          // barycentric interpolation of normal
          const auto& normals = mesh_->normals();
          const auto& normal_indices = mesh_->normal_indices();
          shading_normal_w = w0 * normals[normal_indices[fid][0]] +
                             w1 * normals[normal_indices[fid][1]] +
                             w2 * normals[normal_indices[fid][2]];
        }

        // set shading normal
        if (normal != nullptr) {
          Eigen::Vector3f shading_normal_c =
              w2c_R * shading_normal_w;  // rotate to camera coordinate
          for (int k = 0; k < 3; k++) {
            normal->at(x, y, k) = shading_normal_c[k];
          }
        }

        // delegate color calculation to pixel_shader
        if (color != nullptr) {
          Eigen::Vector3f light_dir = ray_w;  // emit light same as ray
          PixelShaderInput pixel_shader_input(color, x, y, w1, w2, fid, &ray_w,
                                              &light_dir, &shading_normal_w,
                                              &oren_nayar_param, mesh_);
          pixel_shader->Process(pixel_shader_input);
        }
      }
    }
  }

  timer.End();
  LOGI("  Rendering main loop time: %.1f msecs\n", timer.elapsed_msec());

  return true;
}

bool Rasterizer::Impl::RenderColor(Image3b* color) const {
  return Render(color, nullptr, nullptr, nullptr, nullptr);
}

bool Rasterizer::Impl::RenderDepth(Image1f* depth) const {
  return Render(nullptr, depth, nullptr, nullptr, nullptr);
}

bool Rasterizer::Impl::RenderNormal(Image3f* normal) const {
  return Render(nullptr, nullptr, normal, nullptr, nullptr);
}

bool Rasterizer::Impl::RenderMask(Image1b* mask) const {
  return Render(nullptr, nullptr, nullptr, mask, nullptr);
}

bool Rasterizer::Impl::RenderFaceId(Image1i* face_id) const {
  return Render(nullptr, nullptr, nullptr, nullptr, face_id);
}

bool Rasterizer::Impl::RenderW(Image3b* color, Image1w* depth, Image3f* normal,
                               Image1b* mask, Image1i* face_id) const {
  if (depth == nullptr) {
    LOGE("depth is nullptr");
    return false;
  }

  Image1f f_depth;
  bool org_ret = Render(color, &f_depth, normal, mask, face_id);

  if (org_ret) {
    f_depth.ConvertTo(depth);
  }

  return org_ret;
}

bool Rasterizer::Impl::RenderDepthW(Image1w* depth) const {
  return RenderW(nullptr, depth, nullptr, nullptr, nullptr);
}

// Renderer implementation
Rasterizer::Rasterizer() : pimpl_(std::unique_ptr<Impl>(new Impl)) {}

Rasterizer::~Rasterizer() {}

Rasterizer::Rasterizer(const RendererOption& option)
    : pimpl_(std::unique_ptr<Impl>(new Impl(option))) {}

void Rasterizer::set_option(const RendererOption& option) {
  pimpl_->set_option(option);
}

void Rasterizer::set_mesh(std::shared_ptr<const Mesh> mesh) {
  pimpl_->set_mesh(mesh);
}

bool Rasterizer::PrepareMesh() { return pimpl_->PrepareMesh(); }

void Rasterizer::set_camera(std::shared_ptr<const Camera> camera) {
  pimpl_->set_camera(camera);
}

bool Rasterizer::Render(Image3b* color, Image1f* depth, Image3f* normal,
                        Image1b* mask, Image1i* face_id) const {
  return pimpl_->Render(color, depth, normal, mask, face_id);
}

bool Rasterizer::RenderColor(Image3b* color) const {
  return pimpl_->RenderColor(color);
}

bool Rasterizer::RenderDepth(Image1f* depth) const {
  return pimpl_->RenderDepth(depth);
}

bool Rasterizer::RenderNormal(Image3f* normal) const {
  return pimpl_->RenderNormal(normal);
}

bool Rasterizer::RenderMask(Image1b* mask) const {
  return pimpl_->RenderMask(mask);
}

bool Rasterizer::RenderFaceId(Image1i* face_id) const {
  return pimpl_->RenderFaceId(face_id);
}

bool Rasterizer::RenderW(Image3b* color, Image1w* depth, Image3f* normal,
                         Image1b* mask, Image1i* face_id) const {
  return pimpl_->RenderW(color, depth, normal, mask, face_id);
}

bool Rasterizer::RenderDepthW(Image1w* depth) const {
  return pimpl_->RenderDepthW(depth);
}

}  // namespace currender
