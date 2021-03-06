/*
 * Copyright (C) 2019, unclearness
 * All rights reserved.
 */

#include <stdio.h>
#include <fstream>
#include <iostream>
#include <vector>

#include "currender/rasterizer.h"
#include "currender/raytracer.h"
#include "ugu/util.h"

// #define USE_RASTERIZER

using currender::Camera;
using currender::Depth2Gray;
using currender::Depth2Mesh;
using currender::Depth2PointCloud;
using currender::FaceId2RandomColor;
using currender::Image1b;
using currender::Image1f;
using currender::Image1i;
using currender::Image1w;
using currender::Image3b;
using currender::Image3f;
using currender::imwrite;
using currender::Mesh;
using currender::MeshStats;
using currender::Normal2Color;
using currender::PinholeCamera;
using currender::Renderer;
using currender::RendererOption;
using currender::WriteFaceIdAsText;

namespace {
void Test(const std::string& out_dir, std::shared_ptr<Mesh> mesh,
          std::shared_ptr<Camera> camera, const Renderer& renderer) {
  // images
  Image3b color;
  Image1f depth;
  Image1w depthw;
  Image3f normal;
  Image1b mask;
  Image1i face_id;
  Image1b vis_depth;
  Image3b vis_normal;
  Image3b vis_face_id;
  Mesh view_mesh, view_point_cloud;
  const float kMaxConnectZDiff = 100.0f;

  // for pose output by tum format
  std::vector<Eigen::Affine3d> poses;

  MeshStats stats = mesh->stats();
  Eigen::Vector3f center = stats.center;
  Eigen::Vector3f eye;
  Eigen::Matrix4f c2w_mat;

  // translation offset is the largest edge length of bounding box * 1.5
  Eigen::Vector3f diff = stats.bb_max - stats.bb_min;
  float offset = std::max(diff[0], std::max(diff[1], diff[2])) * 1.5f;

  std::vector<Eigen::Affine3d> pose_list;
  std::vector<std::string> name_list;

  // from front
  eye = center;
  eye[2] -= offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, -1, 0), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("front");

  // from back
  eye = center;
  eye[2] += offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, -1, 0), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("back");

  // from right
  eye = center;
  eye[0] += offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, -1, 0), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("right");

  // from left
  eye = center;
  eye[0] -= offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, -1, 0), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("left");

  // from top
  eye = center;
  eye[1] -= offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, 0, 1), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("top");

  // from bottom
  eye = center;
  eye[1] += offset;
  currender::c2w(eye, center, Eigen::Vector3f(0, 0, -1), &c2w_mat);
  pose_list.push_back(Eigen::Affine3d(c2w_mat.cast<double>()));
  name_list.push_back("bottom");

  for (size_t i = 0; i < pose_list.size(); i++) {
    Eigen::Affine3d& c2w = pose_list[i];
    std::string& prefix = name_list[i];

    camera->set_c2w(c2w);
    renderer.Render(&color, &depth, &normal, &mask, &face_id);
    imwrite(out_dir + prefix + "_color.png", color);
    imwrite(out_dir + prefix + "_mask.png", mask);
    currender::ConvertTo(depth, &depthw);
    imwrite(out_dir + prefix + "_depth.png", depthw);
    Depth2Gray(depth, &vis_depth);
    imwrite(out_dir + prefix + "_vis_depth.png", vis_depth);
    Normal2Color(normal, &vis_normal);
    imwrite(out_dir + prefix + "_vis_normal.png", vis_normal);
    FaceId2RandomColor(face_id, &vis_face_id);
    imwrite(out_dir + prefix + "_vis_face_id.png", vis_face_id);
    WriteFaceIdAsText(face_id, out_dir + prefix + "_face_id.txt");
    Depth2Mesh(depth, color, *camera, &view_mesh, kMaxConnectZDiff);
    Depth2PointCloud(depth, color, *camera, &view_point_cloud);
    view_point_cloud.WritePly(out_dir + prefix + "_mesh.ply");
    view_mesh.WriteObj(out_dir, prefix + "_mesh");
    poses.push_back(camera->c2w());
  }

  currender::WriteTumFormat(poses, out_dir + "tumpose.txt");
}

void AlignMesh(std::shared_ptr<Mesh> mesh) {
  // move center as origin
  MeshStats stats = mesh->stats();
  mesh->Translate(-stats.center);

  // rotate 180 deg. around x_axis to align z:forward, y:down, x:right,
  Eigen::Matrix3f R =
      Eigen::AngleAxisf(currender::radians(180.0f), Eigen::Vector3f::UnitX())
          .toRotationMatrix();
  mesh->Rotate(R);

  // recover original translation
  mesh->Translate(stats.center);
}
}  // namespace

int main(int argc, char* argv[]) {
  (void)argc;
  (void)argv;

  // CMake downloads and unzip this data automatically
  // Please do manually if it got wrong
  // http://www.kunzhou.net/tex-models/bunny.zip
  /*
   *  @article{Texturemontage05,
   *  author = "Kun Zhou and Xi Wang and Yiying Tong and Mathieu Desbrun and
   *  Baining Guo and Heung-Yeung Shum", title = "Texturemontage: Seamless
   *  Texturing of Arbitrary Surfaces From Multiple Images", journal = "ACM
   *  Transactions on Graphics", volume = "24", number = "3", year="2005", pages
   *  = "1148-1155"
   */
  std::string data_dir = "../data/bunny/";
  std::string obj_path = data_dir + "bunny.obj";

  std::ifstream ifs(obj_path);
  if (!ifs.is_open()) {
    printf("Please put %s\n", obj_path.c_str());
    return -1;
  }

  // load mesh
  std::shared_ptr<Mesh> mesh = std::make_shared<Mesh>();
  mesh->LoadObj(obj_path, data_dir);

  // original mesh with z:backward, y:up, x:right, like OpenGL
  // align z:forward, y:down, x:right
  AlignMesh(mesh);

  // initialize renderer with diffuse texture color and lambertian shading
  RendererOption option;
  option.diffuse_color = currender::DiffuseColor::kTexture;
  option.diffuse_shading = currender::DiffuseShading::kLambertian;
#ifdef USE_RASTERIZER
  std::unique_ptr<currender::Renderer> renderer =
      std::make_unique<currender::Rasterizer>(option);
#else
  std::unique_ptr<currender::Renderer> renderer =
      std::make_unique<currender::Raytracer>(option);
#endif

  // set mesh
  renderer->set_mesh(mesh);

  // prepare mesh for rendering (e.g. make BVH)
  renderer->PrepareMesh();

  // Make PinholeCamera
  // borrow KinectV1 intrinsics of Freiburg 1 RGB
  // https://vision.in.tum.de/data/datasets/rgbd-dataset/file_formats
  float r = 0.5f;  // scale to smaller size from VGA
  int width = static_cast<int>(640 * r);
  int height = static_cast<int>(480 * r);
  Eigen::Vector2f principal_point(318.6f * r, 255.3f * r);
  Eigen::Vector2f focal_length(517.3f * r, 516.5f * r);
  std::shared_ptr<Camera> camera = std::make_shared<PinholeCamera>(
      width, height, Eigen::Affine3d::Identity(), principal_point,
      focal_length);

  // set camera
  renderer->set_camera(camera);

  // test
  Test(data_dir, mesh, camera, *renderer);

  return 0;
}
