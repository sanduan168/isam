/**
 * @file Loader.cpp
 * @brief Loading files with constraints/factors.
 * @author Michael Kaess
 * @version $Id: Loader.cpp 6902 2012-06-26 02:43:17Z kaess $
 *
 * Copyright (C) 2009-2013 Massachusetts Institute of Technology.
 * Michael Kaess, Hordur Johannsson, David Rosen,
 * Nicholas Carlevaris-Bianco and John. J. Leonard
 *
 * This file is part of iSAM.
 *
 * iSAM is free software; you can redistribute it and/or modify it under
 * the terms of the GNU Lesser General Public License as published by the
 * Free Software Foundation; either version 2.1 of the License, or (at
 * your option) any later version.
 *
 * iSAM is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
 * License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with iSAM.  If not, see <http://www.gnu.org/licenses/>.
 *
 */

#include <iostream>
#include <vector>
#include <map>
#include <list>
#include <utility>

#include <isam/slam2d.h>

#include "Loader.h"

// set to true for external data files with other convention
const bool Y_FORWARD = false;

using namespace std;
using namespace isam;
using namespace Eigen;

/**
 * Create first node at origin: we add a prior to keep the
 * first pose at the origin, which is an arbitrary choice.
 */
void Loader::add_prior() {
  _nodes.resize(1);
  _factors.resize(1);
  _pose_mapper.add(0);
  // add 2D or 3D prior
  if (_is_3d) {
    // create first node
    Pose3d pose0;
    Noise noise = SqrtInformation(100. * eye(6));
    Pose3d_Node* new_pose_node = new Pose3d_Node();
    _nodes[_pose_mapper[0]].push_back(new_pose_node);
    if (_verbose) cout << *new_pose_node << endl;
    _pose_nodes.push_back(new_pose_node);
    // create prior measurement
    Pose3d_Factor* prior = new Pose3d_Factor(new_pose_node, pose0, noise);
    _factors[0].push_back(prior);
    if (_verbose) cout << *prior << endl;
  } else {
    Pose2d pose0;
    Noise noise = SqrtInformation(100. * eye(3));
    Pose2d_Node* new_pose_node = new Pose2d_Node();
    _nodes[_pose_mapper[0]].push_back(new_pose_node);
    if (_verbose) cout << *new_pose_node << endl;
    _pose_nodes.push_back(new_pose_node);
    Pose2d_Factor* prior = new Pose2d_Factor(new_pose_node, pose0, noise);
    _factors[0].push_back(prior);
    if (_verbose) cout << *prior << endl;
  }
}

bool Loader::advance(unsigned int idx_x1, unsigned int next_point_id) {
  // advance to next time step if needed
  bool added = _pose_mapper.add(idx_x1);
  if (added ) {
    _step++;
    _nodes.resize(_step+1);
    _factors.resize(_step+1);
    _num_points.resize(_step+1);
    _num_points[_step] = next_point_id;
    _num_constraints.resize(_step+1);
    _num_constraints[_step] = _constraints.size();
    _num_measurements.resize(_step+1);
    _num_measurements[_step] = _measurements.size();
  }
  return added;
}

void Loader::add_odometry(unsigned int idx_x0, unsigned int idx_x1, const Pose2d& measurement, const Noise& noise) {
  if (advance(idx_x1, _point_nodes.size())) {
    Pose2d_Node* new_pose_node = new Pose2d_Node();
    _nodes[_step].push_back(new_pose_node);
    _pose_nodes.push_back(new_pose_node);
    if (_verbose) cout << idx_x1 << " " << *new_pose_node << endl;
  }
  unsigned int i_x0 = _pose_mapper[idx_x0];
  unsigned int i_x1 = _pose_mapper[idx_x1];
  Pose2d_Pose2d_Factor* factor = new Pose2d_Pose2d_Factor(
      dynamic_cast<Pose2d_Node*>(_pose_nodes[i_x0]),
      dynamic_cast<Pose2d_Node*>(_pose_nodes[i_x1]),
      measurement, noise);
  _factors[i_x1].push_back(factor);
  _constraints.push_back(make_pair(i_x0, i_x1));
  _num_constraints[i_x1] = _constraints.size();
  if (_verbose) cout << i_x1 << " " << *factor << endl;
}

void Loader::add_odometry3(unsigned int idx_x0, unsigned int idx_x1, const Pose3d& measurement, const Noise& noise) {
  if (advance(idx_x1, _point_nodes.size())) {
    Pose3d_Node* new_pose_node = new Pose3d_Node();
    _nodes[_step].push_back(new_pose_node);
    _pose_nodes.push_back(new_pose_node);
    if (_verbose) cout << idx_x1 << " " << *new_pose_node << endl;
  }
  unsigned int i_x0 = _pose_mapper[idx_x0];
  unsigned int i_x1 = _pose_mapper[idx_x1];
  Pose3d_Pose3d_Factor* factor = new Pose3d_Pose3d_Factor(
      dynamic_cast<Pose3d_Node*>(_pose_nodes[i_x0]),
      dynamic_cast<Pose3d_Node*>(_pose_nodes[i_x1]),
      measurement, noise);
  _factors[i_x1].push_back(factor);
  _constraints.push_back(make_pair(i_x0, i_x1));
  _num_constraints[i_x1] = _constraints.size();
  if (_verbose) cout << i_x1 << " " << *factor << endl;
}

void Loader::add_measurement(unsigned int idx_x, unsigned int idx_l, const Point2d& measurement, const Noise& noise) {
  if (_point_mapper.add(idx_l)) {
    // new point has to be added
    Point2d_Node* new_point_node = new Point2d_Node();
    _nodes[_step].push_back(new_point_node);
    _point_nodes.push_back(new_point_node);
    _num_points[_step] = _point_nodes.size();
    if (_verbose) cout << idx_x << " " << *new_point_node << endl;
  }
  unsigned int i_x = _pose_mapper[idx_x];
  unsigned int i_l = _point_mapper[idx_l];
  Pose2d_Point2d_Factor* factor = new Pose2d_Point2d_Factor(
      dynamic_cast<Pose2d_Node*>(_pose_nodes[i_x]),
      dynamic_cast<Point2d_Node*>(_point_nodes[i_l]),
      measurement, noise);
  _factors[i_x].push_back(factor);
  _measurements.push_back(make_pair(i_x, i_l));
  _num_measurements[i_x] = _measurements.size();
  if (_verbose) cout << i_x << " " << *factor << endl;

  for (int i = 0; i < points_landmark.size(); i++) {
//    std::cout << "idx_l: " << idx_l << " " << points_landmark[i].idx << std::endl;
    if (points_landmark[i].idx == idx_l && !points_landmark_prior_flag.at(idx_l)) {
//    if (points_landmark[i].idx == idx_l) {
      Point2d point_2d(points_landmark[i].x,points_landmark[i].y);
      Noise Qsqinf = SqrtInformation(100000000000 * eye(2));
      //Pose2d z0 = x0;
      //z0 = add_noise(z0, sigma);
      //Pose2d_Node *n0 = new Pose2d_Node();
      //Pose2d_Factor *z0f = new Pose2d_Factor(n0, z0, Qsqinf);
      Point2d_Factor *prior = new Point2d_Factor(dynamic_cast<Point2d_Node *>(_point_nodes[i_l]), point_2d, Qsqinf);
      _factors[i_x].push_back(prior);
      points_landmark_prior_flag.at(idx_l) = true;
      break;
    }
  }

}

void Loader::add_point3(unsigned int idx_x, unsigned int idx_p, const Point3d& m, const Noise& noise) {
  // We require that the pose node already exists
  if (_point_mapper.add(idx_p)) {
    Point3d_Node* new_point_node = new Point3d_Node();
    _nodes[_step].push_back(new_point_node);
    _point_nodes.push_back(new_point_node);
    _num_points[_step] = _point_nodes.size();
    if (_verbose) cout << idx_x << " " << *new_point_node << endl;
  }
  unsigned int i_x = _pose_mapper[idx_x];
  unsigned int i_p = _point_mapper[idx_p];
  Pose3d_Point3d_Factor* factor =
    new Pose3d_Point3d_Factor(dynamic_cast<Pose3d_Node*>(_pose_nodes[i_x]),
                              dynamic_cast<Point3d_Node*>(_point_nodes[i_p]),
                              m, noise);
  _factors[i_x].push_back(factor);
  _measurements.push_back(make_pair(i_x, i_p));
  _num_measurements[i_x] = _measurements.size();
  if (_verbose) cout << i_x << " " << *factor << endl;
}

void Loader::add_stereo(isam::StereoCamera* camera,
                        unsigned int idx_x,
                        unsigned int idx_p,
                        const StereoMeasurement& m,
                        const Noise& noise)
{
  // We require that the pose node already exists
  if (_point_mapper.add(idx_p)) {
    Point3dh_Node* new_point_node = new Point3dh_Node();
//    Point3d_Node* new_point_node = new Point3d_Node();
    _nodes[_step].push_back(new_point_node);
    _point_nodes.push_back(new_point_node);
    _num_points[_step] = _point_nodes.size();
    if (_verbose) cout << idx_x << " " << *new_point_node << endl; 
  }
  unsigned int i_x = _pose_mapper[idx_x];
  unsigned int i_p = _point_mapper[idx_p];
  Stereo_Factor* factor = 
    new Stereo_Factor(dynamic_cast<Pose3d_Node*>(_pose_nodes[i_x]),
                      dynamic_cast<Point3dh_Node*>(_point_nodes[i_p]),
//                      dynamic_cast<Point3d_Node*>(_point_nodes[i_p]),
                      camera, m, noise);
  _factors[i_x].push_back(factor);
  _measurements.push_back(make_pair(i_x, i_p));
  _num_measurements[i_x] = _measurements.size();
  if (_verbose) cout << i_x << " " << *factor << endl;
}

bool Loader::parse_line(char* str) {
  bool solve = false;
  char keyword_c[2000];
  int key_length;
  sscanf(str, "%s%n", keyword_c, &key_length);
  const char* arguments = &str[key_length];
  string keyword(keyword_c);
  if (keyword == "ODOMETRY" || keyword == "EDGE2" || keyword == "Constraint2") {
    unsigned int idx_x0, idx_x1;
    double x, y, t, ixx, ixy, ixt, iyy, iyt, itt;
    int res = sscanf(arguments, "%i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg", &idx_x0, &idx_x1, &x, &y, &t, &ixx, &ixy, &ixt, &iyy, &iyt, &itt);
    if (res!=11) {
      cout << "Error while parsing ODOMETRY entry" << endl;
      exit(1);
    }
#if 0 // generate ground truth
    x = round(x);
    y = round(y);
    t = round(t/(M_PI/2)) * (M_PI/2);
    printf("EDGE2 %i %i %g %g %g 50 0 0 50 0 100\n", idx_x0, idx_x1, x, y, t);
#endif
    Pose2d measurement(x, y, t);
    MatrixXd sqrtinf(3,3);
    sqrtinf <<
      ixx, ixy, ixt,
      0.,  iyy, iyt,
      0.,   0., itt;
    if (_step==0) {
      add_prior();
    }
    add_odometry(idx_x0, idx_x1, measurement, SqrtInformation(sqrtinf));
  } else if (keyword == "LANDMARK") {
    unsigned int idx_x, idx_l;
    double x, y, ixx, ixy, iyy;
    int res = sscanf(arguments, "%i %i %lg %lg %lg %lg %lg", &idx_x, &idx_l, &x, &y, &ixx, &ixy, &iyy);
    if (res!=7) {
      cout << "Error while parsing LANDMARK entry" << endl;
      exit(1);
    }
#if 0 // generate ground truth
    x = round(x+0.5)-0.5;
    y = round(y+0.5)-0.5;
    printf("LANDMARK %i %i %g %g 10 0 10\n", idx_x, idx_l, x, y);
#endif
    Point2d measurement(x, y);
    MatrixXd sqrtinf(2,2);
    sqrtinf <<
      ixx, ixy,
      0.,  iyy;
    add_measurement(idx_x, idx_l, measurement, SqrtInformation(sqrtinf));
  } else if (keyword == "POINT3") {
    unsigned int idx_x; // camera node
    unsigned int idx_p; // point node
    double x,y,z; // Measurement, Euclidean coordinates in camera frame
    double i11, i12, i13, i22, i23, i33;
    int res = sscanf(arguments, "%i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg",&idx_x,&idx_p,&x,&y,&z,&i11,&i12,&i13,&i22,&i23,&i33);

    MatrixXd sqrtinf(3,3);
    if (res!=11 && res!=5) {
      cout << "Error while parsing POINT3 entry" << endl;
      exit(1);
    }
    if (res == 5) {
      sqrtinf = eye(3); // no information matrix: use identity
    } else {
      sqrtinf <<
       i11, i12, i13,
         0, i22, i23,
         0,   0, i33;
    }
    Point3d m(x,y,z);
    add_point3(idx_x, idx_p, m, SqrtInformation(sqrtinf));
  } else if (keyword == "CAMERA_STEREO") {
    unsigned int idx_camera; // camera id
    double f;
    double cx, cy;
    double b;

    int res = sscanf(arguments, "%i %lg %lg %lg %lg",&idx_camera,&f,&cx,&cy,&b);
    if (res != 5) {
      cout << "Error while parsing CAMERA_STEREO entry" << endl;
      exit(1);
    }

    _cameras[idx_camera] = new StereoCamera(f, Vector2d(cx,cy), b);

  } else if (keyword == "POINT_STEREO") {
    unsigned int idx_camera; // camera id
    unsigned int idx_x; // camera node
    unsigned int idx_p; // point node
    double u,v,u2; // Measurement, where w is disparity
    double i11, i12, i13, i22, i23, i33;
    int res = sscanf(arguments, "%i %i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg",&idx_camera,&idx_p,&idx_x,&u,&v,&u2,&i11,&i12,&i13,&i22,&i23,&i33);

    MatrixXd sqrtinf(3,3);
    if (res!=12 && res!=6) {
      cout << "Error while parsing POINT_STEREO entry" << endl;
      exit(1);
    }
    if (res == 6) {
      sqrtinf = 4.0 * eye(3); // no information matrix: use identity
    } else {
      sqrtinf <<
       i11, i12, i13,
         0, i22, i23,
         0,   0, i33;
    }
    StereoMeasurement m(u,v,u2);
    if (_cameras.find(idx_camera) == _cameras.end()) {
      cout << "Error while parsing POINT_STEREO entry: Camera (" << idx_camera << ") calibration not found." << endl;
      exit(1);
    }
    add_stereo(_cameras[idx_camera], idx_x, idx_p, m, SqrtInformation(sqrtinf));
  } else if (keyword == "EDGE3") {
    unsigned int idx_x0, idx_x1;
    double x, y, z, yaw, pitch, roll, i11, i12, i13, i14, i15, i16;
    double i22, i23, i24, i25, i26, i33, i34, i35, i36, i44, i45, i46, i55, i56, i66;
    int res = sscanf(arguments, "%i %i %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg %lg",
                     &idx_x0, &idx_x1, &x, &y, &z, &roll, &pitch, &yaw, // note reverse order of angles, also see covariance below
                     &i11, &i12, &i13, &i14, &i15, &i16, &i22, &i23, &i24, &i25, &i26,
                     &i33, &i34, &i35, &i36, &i44, &i45, &i46, &i55, &i56, &i66);
    if (res!=29 && res!=8) {
      cout << "Error while parsing EDGE3 entry" << endl;
      exit(1);
    }
    Pose3d delta;
    if (Y_FORWARD) {
      delta = Pose3d(y, -x, z, yaw, pitch, roll); // converting from external format with Y pointing forward
    } else {
      delta = Pose3d(x, y, z, yaw, pitch, roll); // standard convention (X points forward)
    }
    // square root information matrix
    MatrixXd sqrtinf(6,6);
    if (res==8) {
      sqrtinf = eye(6); // no information matrix: use identity
    } else {
      sqrtinf <<
        i11, i12, i13, i14, i15, i16,
        0., i22, i23, i24, i25, i26,
        0.,  0., i33, i34, i35, i36,
        //                0.,  0.,  0., i44, i45, i46,
        //                0.,  0.,  0.,  0., i55, i56,
        //                0.,  0.,  0.,  0.,  0., i66);
        // todo: this is wrong in the presence of off-diagonal entries because sqrtinf...
        0.,  0.,  0., i66, i56, i46, // note: reversed yaw, pitch, roll, also see sscanf above
        0.,  0.,  0.,  0., i55, i45,
        0.,  0.,  0.,  0.,  0., i44;
    }
    // reverse constraint if needed
    unsigned int i, j;
    if (idx_x0<idx_x1) {
      i = idx_x1;
      j = idx_x0;
    } else {
      delta = Pose3d(delta.oTw());
      i = idx_x0;
      j = idx_x1;
    }
    if (_step==0) {
      _is_3d = true;
      add_prior();
    }
    add_odometry3(j, i, delta, SqrtInformation(sqrtinf));
  } else if (keyword == "POSE3D_INIT") {
    // Initialize a POSE3D node
    unsigned int idx_x0;
    double x, y, z, yaw, pitch, roll;
    int res = sscanf(arguments, "%i %lg %lg %lg %lg %lg %lg",
                     &idx_x0, &x, &y, &z, &roll, &pitch, &yaw); // note reverse order of angles, also see covariance below
    if (res!=7) {
      cout << "Error while parsing POSE3D_INIT entry" << endl;
      exit(1);
    }

    if (advance(idx_x0, _point_nodes.size())) {
      Pose3d_Node* new_pose_node = new Pose3d_Node();
      _nodes[_step].push_back(new_pose_node);
      _pose_nodes.push_back(new_pose_node);
      if (_verbose) cout << idx_x0 << " " << *new_pose_node << endl;
      new_pose_node->init(Pose3d(x,y,z,yaw,pitch,roll));
    }

  } else if (keyword == "EDGE3_INIT") {
    // Provide an edge that is only used for initialization.
  } else if (keyword == "POSE3D_TRUE") {
    // Use to calculate a performance metric

  } else if (keyword == "EDGE3_TRUE") {
    // Use to calculate a performance metric

  } else if (keyword == "SOLVE") {
    solve = true;
  }
  return solve;
}

Loader::Loader(const char* fname, int num_lines, bool verbose) {
  _verbose = verbose;
  _step = 0;
  _is_3d = false;

  // parse and process data file
  _in = fopen(fname, "r");
  if (!_in) {
    printf("ERROR: Failed to open log file %s.\n", fname);
    exit(1);
  }
  int i = 0;
  //load landmark prior
//#if 1
//  points_landmark.push_back({5 ,11.6348, -3.20221});
//  points_landmark.push_back({9 ,15.7621, 4.67216});
//  points_landmark.push_back({32 ,27.8734, 6.01341});
//  points_landmark.push_back({34 ,27.7468, -2.056});
//  points_landmark.push_back({41 ,31.2097, -8.69285});
//  points_landmark.push_back({75 ,40.64, -3.89714});
//  points_landmark.push_back({92 ,45.2068, 5.34749});
//  points_landmark.push_back({97 ,45.6973, 11.6384});
//  points_landmark.push_back({103 ,55.3615, -4.13706});
//  points_landmark.push_back({108 ,58.4193, 10.4297});
//  points_landmark.push_back({114 ,62.9817, -9.69014});
//  points_landmark.push_back({117 ,50.7083, -17.8858});
//  points_landmark.push_back({135 ,45.7325, -24.1333});
//  points_landmark.push_back({140 ,56.4939, -33.939});
//  points_landmark.push_back({151 ,47.2408, -31.8643});
//  points_landmark.push_back({158 ,39.535, -13.4676});
//  points_landmark.push_back({164 ,32.8259, -32.1803});
//  points_landmark.push_back({175 ,27.2307, -18.0543});
//  points_landmark.push_back({179 ,29.2392, -10.8349});
//  points_landmark.push_back({189 ,25.6136, -4.26856});
//  points_landmark.push_back({200 ,17.8183, -15.2462});
//  points_landmark.push_back({212 ,2.99363, -13.7295});
//  points_landmark.push_back({249 ,-5.34256, -6.10593});
//  points_landmark.push_back({259 ,-10.1878, -27.8644});
//  points_landmark.push_back({289 ,-22.5408, -26.8577});
//  points_landmark.push_back({307 ,-25.3481, -18.2685});
//  points_landmark.push_back({316 ,-20.5612, -37.1532});
//  points_landmark.push_back({318 ,-18.0613, -37.4471});
//  points_landmark.push_back({320 ,-15.3035, -37.8897});
//  points_landmark.push_back({346 ,-23.4147, -4.57647});
//  points_landmark.push_back({355 ,-34.33, -25.2644});
//  points_landmark.push_back({358 ,-33.0591, -4.42294});
//  points_landmark.push_back({379 ,-20.1737, 3.3986});
//  points_landmark.push_back({383 ,-26.681, 3.12312});
//  points_landmark.push_back({408 ,-42.2504, -2.41835});
//  points_landmark.push_back({436 ,-59.6066, 0.920812});
//  points_landmark.push_back({451 ,-61.4558, -5.81254});
//  points_landmark.push_back({455 ,-41.6377, -10.9096});
//  points_landmark.push_back({466 ,-51.6986, -21.6325});
//  points_landmark.push_back({584 ,48.9734, -25.745});
//  points_landmark.push_back({609 ,46.3556, -55.8904});
//  points_landmark.push_back({632 ,42.6894, -72.05});
//  points_landmark.push_back({636 ,53.2911, -67.6578});
//  points_landmark.push_back({639 ,40.2628, -78.1715});
//  points_landmark.push_back({646 ,48.0989, -71.5668});
//  points_landmark.push_back({661 ,55.9158, -65.2969});
//  points_landmark.push_back({705 ,64.7284, -45.1483});
//  points_landmark.push_back({707 ,82.641, -52.4499});
//  points_landmark.push_back({711 ,72.8516, -38.6782});
//  points_landmark.push_back({724 ,77.8836, -41.7807});
//  points_landmark.push_back({756 ,62.7556, -26.7799});
//  points_landmark.push_back({907 ,2.51708, -39.7508});
//  points_landmark.push_back({1033 ,64.3075, -80.9979});
//  points_landmark.push_back({1039 ,56.7892, -93.6778});
//  points_landmark.push_back({1043 ,50.6411, -97.2419});
//  points_landmark.push_back({1058 ,51.4938, -107.95});
//  points_landmark.push_back({1060 ,37.5494, -103.289});
//  points_landmark.push_back({1071 ,30.8341, -104.531});
//  points_landmark.push_back({1079 ,36.1725, -108.118});
//  points_landmark.push_back({1095 ,54.1603, -124.746});
//  points_landmark.push_back({1112 ,35.0084, -143.702});
//  points_landmark.push_back({1132 ,37.8604, -155.005});
//  points_landmark.push_back({1147 ,40.9888, -166.749});
//  points_landmark.push_back({1154 ,49.9827, -161.196});
//  points_landmark.push_back({1163 ,58.3202, -159.864});
//  points_landmark.push_back({1178 ,66.2733, -164.087});
//  points_landmark.push_back({1193 ,69.578, -179.702});
//  points_landmark.push_back({1203 ,79.6443, -174.516});
//  points_landmark.push_back({1235 ,91.1395, -165.927});
//  points_landmark.push_back({1268 ,96.5601, -162.096});
//  points_landmark.push_back({1309 ,68.7071, -146.118});
//  points_landmark.push_back({1329 ,70.4407, -128.669});
//  points_landmark.push_back({1333 ,61.4379, -140.561});
//  points_landmark.push_back({1389 ,40.3709, -110.573});
//  points_landmark.push_back({1867 ,44.3882, -22.105});
//  points_landmark.push_back({1876 ,56.6894, -41.511});
//  points_landmark.push_back({1878 ,42.8754, -25.1366});
//  points_landmark.push_back({2574 ,182.799, -44.8861});
//  points_landmark.push_back({2758 ,84.373, 6.31319});
//  points_landmark.push_back({3254 ,-2.15893, -28.4063});
//  points_landmark.push_back({3527 ,91.2592, -43.535});
//  points_landmark.push_back({3538 ,86.33, -38.6902});
//  points_landmark.push_back({3906 ,7.21193, 81.4774});
//  points_landmark.push_back({3953 ,-6.29528, 85.1181});
//  points_landmark.push_back({3986 ,-18.96, 89.1166});
//  points_landmark.push_back({4035 ,-16.9574, 87.3557});
//  points_landmark.push_back({4085 ,-11.7848, 48.3019});
//  points_landmark.push_back({4580 ,-123.283, -154.826});
//  points_landmark.push_back({4592 ,-133.541, -146.455});
//  points_landmark.push_back({4605 ,-130.679, -141.459});
//  points_landmark.push_back({4607 ,-132.328, -136.965});
//  points_landmark.push_back({4613 ,-141.364, -154.898});
//  points_landmark.push_back({4628 ,-144.99, -135.683});
//  points_landmark.push_back({4670 ,-160.401, -113.827});
//  points_landmark.push_back({4688 ,-151.41, -93.3861});
//  points_landmark.push_back({4703 ,-143.877, -83.161});
//  points_landmark.push_back({4715 ,-132.796, -77.9779});
//  points_landmark.push_back({4743 ,-163.167, -62.404});
//  points_landmark.push_back({4753 ,-152.036, -48.9739});
//  points_landmark.push_back({4762 ,-161.075, -46.9695});
//  points_landmark.push_back({4788 ,-160.059, -31.7529});
//  points_landmark.push_back({4792 ,-156.198, -21.8899});
//  points_landmark.push_back({4795 ,-131.154, -25.4328});
//  points_landmark.push_back({4797 ,-144.155, -24.582});
//  points_landmark.push_back({4820 ,-159.122, -0.556096});
//  points_landmark.push_back({4826 ,-159.751, 4.91789});
//  points_landmark.push_back({4835 ,-154.579, 10.2052});
//  points_landmark.push_back({4838 ,-154.917, 18.0524});
//  points_landmark.push_back({4848 ,-165.356, 16.1931});
//  points_landmark.push_back({4856 ,-150.983, 30.9938});
//  points_landmark.push_back({4868 ,-139.413, 34.0644});
//  points_landmark.push_back({4876 ,-134.374, 31.9307});
//  points_landmark.push_back({4886 ,-126.219, 37.5167});
//  points_landmark.push_back({4900 ,-121.064, 35.5573});
//  points_landmark.push_back({4961 ,-142.05, 1.197});
//  points_landmark.push_back({4969 ,-138.444, -16.26});
//  points_landmark.push_back({4983 ,-142.253, -22.1611});
//  points_landmark.push_back({4985 ,-138.713, -16.0793});
//  points_landmark.push_back({5446 ,-38.9491, -183.381});
//  points_landmark.push_back({5450 ,-47.0894, -187.231});
//  points_landmark.push_back({5624 ,-135.654, -138.853});
//  points_landmark.push_back({5627 ,-142.352, -122.769});
//  points_landmark.push_back({5872 ,-47.5657, 11.5977});
//  points_landmark.push_back({5878 ,-35.7959, 2.29276});
//  points_landmark.push_back({5884 ,-31.979, -2.59896});
//  points_landmark.push_back({5887 ,-42.9256, 13.5451});
//  points_landmark.push_back({5890 ,-48.721, 20.7677});
//  points_landmark.push_back({5896 ,-35.8173, 15.6358});
//  points_landmark.push_back({5898 ,-29.6545, 11.1612});
//  points_landmark.push_back({5903 ,-31.3334, 2.53848});
//  points_landmark.push_back({5913 ,-38.5297, 21.6964});
//  points_landmark.push_back({5916 ,-24.0601, 16.1066});
//  points_landmark.push_back({5927 ,-21.3001, 9.90365});
//  points_landmark.push_back({5941 ,-10.2912, 19.898});
//  points_landmark.push_back({5973 ,5.72065, 23.5951});
//  points_landmark.push_back({5994 ,6.42425, 45.4804});
//  points_landmark.push_back({6003 ,18.6463, 46.7597});
//  points_landmark.push_back({6008 ,10.6991, 46.9778});
//  points_landmark.push_back({6014 ,28.453, 35.6148});
//  points_landmark.push_back({6019 ,25.4801, 44.5675});
//  points_landmark.push_back({6035 ,30.2317, 43.3014});
//  points_landmark.push_back({6039 ,33.6699, 29.422});
//  points_landmark.push_back({6042 ,38.59, 23.7103});
//  points_landmark.push_back({6068 ,43.1052, 17.9344});
//  points_landmark.push_back({6184 ,-20.2737, 7.44344});
//  points_landmark.push_back({6210 ,-40.8915, 21.0223});
//  points_landmark.push_back({6218 ,-36.8372, 28.1965});
//  points_landmark.push_back({6222 ,-38.5735, 23.8847});
//  points_landmark.push_back({6253 ,-57.1639, 21.4694});
//  points_landmark.push_back({6873 ,-174.544, -189.762});
//  points_landmark.push_back({6884 ,-183.35, -197.313});
//#elif
  points_landmark.push_back({5, 11.5465, -3.17907});
  points_landmark.push_back({9, 15.7492, 4.58071});
  points_landmark.push_back({32, 27.8274, 5.8112});
  points_landmark.push_back({34, 27.4151, -2.28382});
  points_landmark.push_back({41, 30.9937, -8.74633});
  points_landmark.push_back({75, 40.5417, -3.7835});
  points_landmark.push_back({92, 45.162, 5.63353});
  points_landmark.push_back({97, 45.4695, 12.161});
  points_landmark.push_back({103, 55.135, -3.36122});
  points_landmark.push_back({108, 58.5388, 11.1018});
  points_landmark.push_back({114, 62.8829, -8.8011});
  points_landmark.push_back({117, 50.7583, -17.747});
  points_landmark.push_back({135, 46.1381, -24.0809});
  points_landmark.push_back({140, 57.9534, -33.2087});
  points_landmark.push_back({151, 48.0571, -31.7115});
  points_landmark.push_back({158, 40.1351, -12.9558});
  points_landmark.push_back({164, 33.3495, -31.461});
  points_landmark.push_back({175, 27.9234, -16.6886});
  points_landmark.push_back({179, 30.6231, -9.48371});
  points_landmark.push_back({189, 27.3647, -2.88626});
  points_landmark.push_back({200, 18.8215, -13.1128});
  points_landmark.push_back({212, 4.31242, -10.887});
  points_landmark.push_back({249, -4.093, -4.21137});
  points_landmark.push_back({259, -9.01129, -25.9345});
  points_landmark.push_back({289, -21.0757, -24.8959});
  points_landmark.push_back({307, -23.8781, -16.4912});
  points_landmark.push_back({316, -18.921, -35.376});
  points_landmark.push_back({318, -16.4172, -35.6403});
  points_landmark.push_back({320, -13.654, -36.0508});
  points_landmark.push_back({346, -22.5566, -3.19301});
  points_landmark.push_back({355, -32.7949, -23.4161});
  points_landmark.push_back({358, -31.626, -2.66164});
  points_landmark.push_back({379, -18.1323, 4.80523});
  points_landmark.push_back({383, -24.6711, 4.98441});
  points_landmark.push_back({408, -40.5742, -1.87544});
  points_landmark.push_back({436, -57.9772, -0.996518});
  points_landmark.push_back({451, -58.5032, -8.37915});
  points_landmark.push_back({455, -38.2485, -10.3609});
  points_landmark.push_back({466, -45.1999, -22.2852});
  points_landmark.push_back({584, 51.9624, 12.1937});
  points_landmark.push_back({609, 61.7042, -16.0292});
  points_landmark.push_back({632, 66.0112, -32.1951});
  points_landmark.push_back({636, 73.6034, -25.005});
  points_landmark.push_back({639, 65.7051, -38.5328});
  points_landmark.push_back({646, 70.803, -28.6906});
  points_landmark.push_back({661, 75.018, -19.8815});
  points_landmark.push_back({705, 72.0482, 2.34355});
  points_landmark.push_back({707, 91.1066, 5.3867});
  points_landmark.push_back({711, 75.3074, 11.8681});
  points_landmark.push_back({724, 81.6599, 12.1995});
  points_landmark.push_back({756, 58.6899, 11.0758});
  points_landmark.push_back({907, 43.0972, -35.2945});
  points_landmark.push_back({1033, 90.7602, 23.6827});
  points_landmark.push_back({1039, 103.169, 16.0744});
  points_landmark.push_back({1043, 106.738, 10.0124});
  points_landmark.push_back({1058, 117.432, 11.0104});
  points_landmark.push_back({1060, 112.757, -2.72523});
  points_landmark.push_back({1071, 114.145, -9.35412});
  points_landmark.push_back({1079, 117.831, -3.99577});
  points_landmark.push_back({1095, 134.407, 14.1532});
  points_landmark.push_back({1112, 153.725, -4.06898});
  points_landmark.push_back({1132, 164.875, -0.128613});
  points_landmark.push_back({1147, 176.453, 3.81027});
  points_landmark.push_back({1154, 170.139, 12.5881});
  points_landmark.push_back({1163, 168.401, 20.6606});
  points_landmark.push_back({1178, 171.839, 28.9396});
  points_landmark.push_back({1193, 186.401, 34.0626});
  points_landmark.push_back({1203, 179.743, 43.4148});
  points_landmark.push_back({1235, 169.23, 53.4939});
  points_landmark.push_back({1268, 163.88, 57.9345});
  points_landmark.push_back({1309, 154.403, 26.3912});
  points_landmark.push_back({1329, 137.553, 23.2422});
  points_landmark.push_back({1333, 151.298, 17.8505});
  points_landmark.push_back({1389, 130.102, -11.4927});
  points_landmark.push_back({1867, 55.58, -37.709});
  points_landmark.push_back({1876, 62.9993, -15.7718});
  points_landmark.push_back({1878, 59.0859, -37.0088});
  points_landmark.push_back({2574, -44.5505, 6.99611});
  points_landmark.push_back({ 2758, 53.2813, -40.2982 });
  points_landmark.push_back({ 3254, 119.628, 19.6877 });
  points_landmark.push_back({ 3527, 71.6752, -24.9713 });
  points_landmark.push_back({ 3538, 76.739, -20.1949 });
  points_landmark.push_back({ 3906, 188.905, 66.1358 });
  points_landmark.push_back({ 3953, 188.681, 80.2309 });
  points_landmark.push_back({ 3986, 188.847, 92.706 });
  points_landmark.push_back({ 4035, 187.495, 88.2711 });
  points_landmark.push_back({ 4085, 156.602, 63.0981 });
  points_landmark.push_back({ 4580, -75.8723, -0.283061 });
  points_landmark.push_back({ 4592, -79.0246, 12.3848 });
  points_landmark.push_back({ 4605, -73.9653, 14.714 });
  points_landmark.push_back({ 4607, -72.571, 19.0486 });
  points_landmark.push_back({ 4613, -90.4048, 9.73471 });
  points_landmark.push_back({ 4628, -82.4706, 27.4213 });
  points_landmark.push_back({ 4670, -82.3299, 54.3869 });
  points_landmark.push_back({ 4688, -63.663, 65.7114 });
  points_landmark.push_back({ 4703, -51.2485, 69.9158 });
  points_landmark.push_back({ 4715, -39.3992, 68.3781 });
  points_landmark.push_back({ 4743, -56.4995, 98.109 });
  points_landmark.push_back({ 4753, -40.2023, 103.439 });
  points_landmark.push_back({ 4762, -47.0425, 109.658 });
  points_landmark.push_back({ 4788, -38.3961, 122.262 });
  points_landmark.push_back({ 4792, -29.9772, 129.971 });
  points_landmark.push_back({ 4795, -10.3363, 113.427 });
  points_landmark.push_back({ 4797, -21.229, 120.636 });
  points_landmark.push_back({ 4820, -22.8838, 148.911 });
  points_landmark.push_back({ 4826, -20.8566, 153.985 });
  points_landmark.push_back({ 4835, -14.4251, 155.52 });
  points_landmark.push_back({ 4838, -10.9666, 163.422 });
  points_landmark.push_back({ 4848, -20.9836, 166.456 });
  points_landmark.push_back({ 4856, -1.78278, 173.201 });
  points_landmark.push_back({ 4868, 9.87043, 171.018 });
  points_landmark.push_back({ 4876, 13.5191, 166.842 });
  points_landmark.push_back({ 4886, 23.3449, 168.39 });
  points_landmark.push_back({ 4900, 27.2061, 164.612 });
  points_landmark.push_back({ 4961, -6.33903, 142.445 });
  points_landmark.push_back({ 4969, -10.2514, 124.986 });
  points_landmark.push_back({ 4983, -17.6765, 119.916 });
  points_landmark.push_back({ 4985, -10.4929, 125.221 });
  points_landmark.push_back({ 5446, 54.4987, -34.518 });
  points_landmark.push_back({ 5450, 45.8073, -38.3984 });
  points_landmark.push_back({ 5624, -43.4568, 9.66857 });
  points_landmark.push_back({ 5627, -51.8568, 24.836 });
  points_landmark.push_back({ 5872, 9.62202, 179.312 });
  points_landmark.push_back({ 5878, 23.7279, 173.974 });
  points_landmark.push_back({ 5884, 28.8615, 170.856 });
  points_landmark.push_back({ 5887, 13.2693, 182.632 });
  points_landmark.push_back({ 5890, 5.46962, 187.518 });
  points_landmark.push_back({ 5896, 19.3689, 186.67 });
  points_landmark.push_back({ 5898, 26.5869, 184.723 });
  points_landmark.push_back({ 5903, 27.7708, 176.108 });
  points_landmark.push_back({ 5913, 14.4592, 191.73 });
  points_landmark.push_back({ 5916, 30.5502, 190.803 });
  points_landmark.push_back({ 5927, 34.8006, 186.666 });
  points_landmark.push_back({ 5941, 41.2741, 200.072 });
  points_landmark.push_back({ 5973, 54.5185, 209.903 });
  points_landmark.push_back({ 5994, 46.1958, 230.511 });
  points_landmark.push_back({ 6003, 56.5828, 236.923 });
  points_landmark.push_back({ 6008, 49.2759, 233.676 });
  points_landmark.push_back({ 6014, 70.1584, 231.199 });
  points_landmark.push_back({ 6019, 63.4056, 237.902 });
  points_landmark.push_back({ 6035, 68.2317, 238.889 });
  points_landmark.push_back({ 6039, 77.4937, 228.03 });
  points_landmark.push_back({ 6042, 84.3808, 225.386 });
  points_landmark.push_back({ 6068, 91.1285, 222.256 });
  points_landmark.push_back({ 6184, 40.1288, 183.843 });
  points_landmark.push_back({ 6210, 15.6948, 185.569 });
  points_landmark.push_back({ 6218, 15.6903, 193.741 });
  points_landmark.push_back({ 6222, 16.3519, 189.08 });
  points_landmark.push_back({ 6253, 1.6775, 177.032 });
  points_landmark.push_back({ 6873, 71.779, -21.6572 });
  points_landmark.push_back({ 6884, 74.7769, -33.0628 });
//#endif
  for (int j = 0; j < points_landmark.size(); j++) {
    points_landmark_prior_flag.insert(make_pair(points_landmark[j].idx,false));
  }
  while (!feof(_in) && (num_lines == 0 || i < num_lines)) {
    char str[2000];
    if (fgets(str, 2000, _in)) {
      parse_line(str);
    }
    i++;
  }
  fclose(_in);
}

Loader::~Loader()
{
  for(std::map<int, isam::StereoCamera*>::iterator it = _cameras.begin(); it != _cameras.end(); ++it) delete it->second;
}

void Loader::print_stats() const {
  int n = num_steps();
  cout << "Number of poses: " << n << endl;
  if (_point_nodes.size()>0) {
    cout << "Number of landmarks: " << _num_points[n-1] << endl;
    cout << "Number of measurements: " << _num_measurements[n-1] << endl;
  }
  cout << "Number of constraints: " << _num_constraints[n-1] << endl;
}

bool Loader::more_data(unsigned int* step) {
  (*step)++;
  return (*step)<=num_steps();
}

const Loader::PoseList Loader::poses(unsigned int step) const {
  Loader::PoseList poses;
  poses.resize(step+1);
  for (unsigned int i=0; i<step+1; i++) {
    if (is_3d()) {
      poses[i] = dynamic_cast<Pose3d_Node*>(_pose_nodes[i])->value();
    } else {
      poses[i].of_pose2d(dynamic_cast<Pose2d_Node*>(_pose_nodes[i])->value());
    }
  }
  return poses;
}

const Loader::PoseList Loader::points(unsigned int step) const {
  Loader::PoseList points;
  points.resize(_num_points[step]);
  for (unsigned int i=0; i<points.size(); i++) {
    if (is_3d()) {
      //points[i].of_point3d(dynamic_cast<Point3d_Node*>(_point_nodes[i])->value());
      if (dynamic_cast<Point3d_Node*>(_point_nodes[i])) {
        points[i].of_point3d(dynamic_cast<Point3d_Node*>(_point_nodes[i])->value());
      } else if (dynamic_cast<Point3dh_Node*>(_point_nodes[i])) {
        Point3dh_Node* node = dynamic_cast<Point3dh_Node*>(_point_nodes[i]);
        points[i].of_point3d(node->value().to_point3d());
      }
    } else {
      points[i].of_point2d(dynamic_cast<Point2d_Node*>(_point_nodes[i])->value());
    }
  }
  return points;
}

const vector<pair<int,int> > Loader::constraints(unsigned int step) const {
  return vector<pair<int,int> >(_constraints.begin(), _constraints.begin()+_num_constraints[step]);
}

const vector<pair<int,int> > Loader::measurements(unsigned int step) const {
  return vector<pair<int,int> >(_measurements.begin(), _measurements.begin()+_num_measurements[step]);
}
