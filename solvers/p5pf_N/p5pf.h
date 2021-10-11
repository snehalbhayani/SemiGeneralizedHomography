/*
 * p35pf.h
 * Author: Viktor Larsson
 */
#pragma once
#include <Eigen/Dense>
#include <vector>

namespace colmap {
	/* Absolute pose with unknown focal length from 5 2D-3D point correspondences */
	int p5pf(const std::vector<Eigen::Vector2d>& points2d,
		const std::vector<Eigen::Vector3d>& points3d,
		std::vector<Eigen::Matrix<double,3,4>>* poses,
		std::vector<double>* focal_lengths,
		bool normalize_input = true);

}