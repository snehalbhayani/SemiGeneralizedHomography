/* MyMEXFunction
 * Adds second input to each
 * element of first input
 * a = MyMEXFunction(a,b);
 */
#pragma once

#include "mex.h"
#include <algorithm>
#include <fstream>
#include <functional>
#include <sstream>
#include <vector>

#include <opencv2/core/core.hpp>
#include <opencv2/highgui/highgui.hpp>
#include <opencv2/imgproc.hpp>
#include <opencv2/features2d.hpp>
#include <Eigen/Eigen>

#include "types.h"
#include "statistics.h"

namespace gcransac
{
	namespace utils
	{
		/*
			Function declaration
		*/
		void drawMatches(
			const cv::Mat &points_,
			const std::vector<size_t> &inliers_,
			const cv::Mat &image1_,
			const cv::Mat &image2_,
			cv::Mat &out_image_,
			int circle_radius_ = 4);

		void drawImagePoints(
			const cv::Mat &points_,
			const std::vector<size_t> &inliers_,
			const cv::Scalar &color_,
			cv::Mat &image_,
			int circle_radius_ = 10);

		bool savePointsToFile(
			const cv::Mat &points_,
			const char* file_,
			const std::vector<size_t> *inliers_ = NULL);

		template<size_t _Dimensions = 4, 
			size_t _StoreEveryKth = 1, 
			bool _PointNumberInHeader = true>
		bool loadPointsFromFile(
			cv::Mat &points_,
			const char* file_);

		void detectFeatures(
			std::string name_,
			cv::Mat image1_,
			cv::Mat image2_,
			cv::Mat &points_);

		void showImage(
			const cv::Mat &image_,
			std::string window_name_,
			size_t max_width_,
			size_t max_height_,
			bool wait_);

		template<typename T, size_t N, size_t M>
		bool loadMatrix(const std::string &path_,
			Eigen::Matrix<T, N, M> &matrix_);

		void normalizeCorrespondences(const cv::Mat &points_,
			const Eigen::Matrix3d &intrinsics_src_,
			const Eigen::Matrix3d &intrinsics_dst_,
			cv::Mat &normalized_points_);

		void saveStatisticsToFile(const gcransac::utils::RANSACStatistics &statistics_,
			const std::string &source_path_,
			const std::string &destination_path_,
			const std::string &filename_,
#ifdef _WIN32
			const int mode_ = std::fstream::app
#else
			const std::ios_base::openmode mode_ = std::fstream::app
#endif
		);

		/*
			Function definition
		*/
		void drawImagePoints(
			const cv::Mat &points_,
			const std::vector<size_t> &inliers_,
			const cv::Scalar &color_,
			cv::Mat &image_,
			int circle_radius_)
		{
			for (const auto &idx : inliers_)
			{
				cv::Point2d point(points_.at<double>(idx, 0),
					points_.at<double>(idx, 1));
				
				cv::circle(image_, 
					point, 
					circle_radius_, 
					color_, 
					static_cast<int>(circle_radius_));
			}
		}

		void drawMatches(
			const cv::Mat &points_,
			const std::vector<size_t> &inliers_,
			const cv::Mat &image_src_,
			const cv::Mat &image_dst_,
			cv::Mat &out_image_,
			int circle_radius_)
		{
			// Final image
			out_image_.create(image_src_.rows, // Height
				2 * image_src_.cols, // Width
				image_src_.type()); // Type

			cv::Mat roi_img_result_left =
				out_image_(cv::Rect(0, 0, image_src_.cols, image_src_.rows)); // Img1 will be on the left part
			cv::Mat roi_img_result_right =
				out_image_(cv::Rect(image_src_.cols, 0, image_dst_.cols, image_dst_.rows)); // Img2 will be on the right part, we shift the roi of img1.cols on the right

			cv::Mat roi_image_src = image_src_;
			cv::Mat roi_image_dst = image_dst_;

			roi_image_src.copyTo(roi_img_result_left); //Img1 will be on the left of imgResult
			roi_image_dst.copyTo(roi_img_result_right); //Img2 will be on the right of imgResult

			for (const auto &idx : inliers_)
			{
				cv::Point2d pt1(points_.at<double>(idx, 0),
					points_.at<double>(idx, 1));
				cv::Point2d pt2(image_dst_.cols + points_.at<double>(idx, 2),
					points_.at<double>(idx, 3));

				cv::Scalar color(255 * static_cast<double>(rand()) / RAND_MAX,
					255 * static_cast<double>(rand()) / RAND_MAX,
					255 * static_cast<double>(rand()) / RAND_MAX);

				cv::circle(out_image_, pt1, circle_radius_, color, static_cast<int>(circle_radius_ * 0.4));
				cv::circle(out_image_, pt2, circle_radius_, color, static_cast<int>(circle_radius_ * 0.4));
				cv::line(out_image_, pt1, pt2, color, 2);
			}
		}

		void detectFeatures(std::string scene_name_,
			cv::Mat image1_,
			cv::Mat image2_,
			cv::Mat &points_)
		{
			if (loadPointsFromFile(points_,
				scene_name_.c_str()))
			{
				printf("Match number: %d\n", points_.rows);
				return;
			}

			printf("Detect AKAZE features\n");
			cv::Mat descriptors1, descriptors2;
			std::vector<cv::KeyPoint> keypoints1, keypoints2;

			cv::Ptr<cv::AKAZE> detector = cv::AKAZE::create();
			detector->detect(image1_, keypoints1);
			detector->compute(image1_, keypoints1, descriptors1);
			printf("Features found in the first image: %d\n", static_cast<int>(keypoints1.size()));

			detector->detect(image2_, keypoints2);
			detector->compute(image2_, keypoints2, descriptors2);
			printf("Features found in the second image: %d\n", static_cast<int>(keypoints2.size()));

			cv::BFMatcher matcher(cv::NORM_HAMMING);
			std::vector< std::vector< cv::DMatch >> matches_vector;
			matcher.knnMatch(descriptors1, descriptors2, matches_vector, 2);

			std::vector<std::tuple<double, cv::Point2d, cv::Point2d>> correspondences;
			for (auto match : matches_vector)
			{
				if (match.size() == 2 && match[0].distance < match[1].distance * 0.8)
				{
					auto& kp1 = keypoints1[match[0].queryIdx];
					auto& kp2 = keypoints2[match[0].trainIdx];
					correspondences.emplace_back(std::make_tuple<double, cv::Point2d, cv::Point2d>(match[0].distance / match[1].distance, (cv::Point2d)kp1.pt, (cv::Point2d)kp2.pt));
				}
			}

			// Sort the points for PROSAC
			std::sort(correspondences.begin(), correspondences.end(), [](const std::tuple<double, cv::Point2d, cv::Point2d>& correspondence_1_,
				const std::tuple<double, cv::Point2d, cv::Point2d>& correspondence_2_) -> bool
			{
				return std::get<0>(correspondence_1_) < std::get<0>(correspondence_2_);
			});

			points_ = cv::Mat(static_cast<int>(correspondences.size()), 4, CV_64F);
			double *points_ptr = reinterpret_cast<double*>(points_.data);

			for (auto[distance_ratio, point_1, point_2] : correspondences)
			{
				*(points_ptr++) = point_1.x;
				*(points_ptr++) = point_1.y;
				*(points_ptr++) = point_2.x;
				*(points_ptr++) = point_2.y;
			}

			savePointsToFile(points_, scene_name_.c_str());
			printf("Match number: %d\n", static_cast<int>(points_.rows));
		}

		template<size_t _Dimensions, size_t _StoreEveryKth, bool _PointNumberInHeader>
		bool loadPointsFromFile(cv::Mat &points,
			const char* file)
		{
			std::ifstream infile(file);

			if (!infile.is_open())
				return false;

			int N;
			std::string line;
			int line_idx = 0, 
				loaded_lines = 0;
			double *points_ptr = NULL;
			cv::Mat point;

			while (getline(infile, line))
			{
				// If the number of points to be loaded is stored in the file's header
				if (_PointNumberInHeader)
				{
					if (line_idx++ == 0)
					{
						N = atoi(line.c_str()) / _StoreEveryKth;
						points = cv::Mat(N, _Dimensions, CV_64F);
						points_ptr = reinterpret_cast<double*>(points.data);
						continue;
					}

					if (line_idx % _StoreEveryKth)
						continue;

					std::istringstream split(line);
					for (size_t dim = 0; dim < _Dimensions; ++dim)
						split >> *(points_ptr++);

					if (++loaded_lines >= N)
						break;
				}
				else // Otherwise
				{
					if (++loaded_lines % _StoreEveryKth)
						continue;
					
					if (point.rows != 1 || point.cols != _Dimensions)
					{
						point.create(1, _Dimensions, CV_64F);
						points_ptr = reinterpret_cast<double*>(point.data);
					}

					std::istringstream split(line);
					for (size_t dim = 0; dim < _Dimensions; ++dim)
						split >> points_ptr[dim];

					// Add the new point as the last row
					if (points.cols != _Dimensions)
						points.create(0, _Dimensions, CV_64F);
					points.push_back(point);
				}
			}

			infile.close();
			return true;
		}

		bool savePointsToFile(const cv::Mat &points,
			const char* file,
			const std::vector<size_t> *inliers)
		{
			std::ofstream outfile(file, std::ios::out);

			double *points_ptr = reinterpret_cast<double*>(points.data);
			const size_t M = points.cols;

			if (inliers == NULL)
			{
				outfile << points.rows << std::endl;
				for (auto i = 0; i < points.rows; ++i)
				{
					for (auto j = 0; j < M; ++j)
						outfile << *(points_ptr++) << " ";
					outfile << std::endl;
				}
			}
			else
			{
				outfile << inliers->size() << std::endl;
				for (size_t i = 0; i < inliers->size(); ++i)
				{
					const size_t offset = inliers->at(i) * M;
					for (auto j = 0; j < M; ++j)
						outfile << *(points_ptr + offset + j) << " ";
					outfile << std::endl;
				}
			}

			outfile.close();

			return true;
		}

		template<typename T, size_t N, size_t M>
		bool loadMatrix(const std::string &path_,
			Eigen::Matrix<T, N, M> &matrix_)
		{
			std::ifstream infile(path_);

			if (!infile.is_open())
				return false;

			size_t row = 0,
				column = 0;
			double element;

			while (infile >> element)
			{
				matrix_(row, column) = element;
				++column;
				if (column >= M)
				{
					column = 0;
					++row;
				}
			}

			infile.close();

			return row == N &&
				column == 0;
		}

		void showImage(const cv::Mat &image_,
			std::string window_name_,
			size_t max_width_,
			size_t max_height_,
			bool wait_)
		{
			// Resizing the window to fit into the screen if needed
			int window_width = image_.cols,
				window_height = image_.rows;
			if (static_cast<double>(image_.cols) / max_width_ > 1.0 &&
				static_cast<double>(image_.cols) / max_width_ >
				static_cast<double>(image_.rows) / max_height_)
			{
				window_width = max_width_;
				window_height = static_cast<int>(window_width * static_cast<double>(image_.rows) / static_cast<double>(image_.cols));
			}
			else if (static_cast<double>(image_.rows) / max_height_ > 1.0 &&
				static_cast<double>(image_.cols) / max_width_ <
				static_cast<double>(image_.rows) / max_height_)
			{
				window_height = max_height_;
				window_width = static_cast<int>(window_height * static_cast<double>(image_.cols) / static_cast<double>(image_.rows));
			}

			cv::namedWindow(window_name_, cv::WINDOW_NORMAL);
			cv::resizeWindow(window_name_, window_width, window_height);
			cv::imshow(window_name_, image_);
			if (wait_)
				cv::waitKey(0);
		}

		void normalizeCorrespondences(const cv::Mat &points_,
			const Eigen::Matrix3d &intrinsics_src_,
			const Eigen::Matrix3d &intrinsics_dst_,
			cv::Mat &normalized_points_)
		{
			const double *points_ptr = reinterpret_cast<double *>(points_.data);
			double *normalized_points_ptr = reinterpret_cast<double *>(normalized_points_.data);
			const Eigen::Matrix3d inverse_intrinsics_src = intrinsics_src_.inverse(),
				inverse_intrinsics_dst = intrinsics_dst_.inverse();

			// Most likely, this is not the fastest solution, but it does
			// not affect the speed of Graph-cut RANSAC, so not a crucial part of
			// this example.
			double x0, y0, x1, y1;
			for (auto r = 0; r < points_.rows; ++r)
			{
				Eigen::Vector3d point_src,
					point_dst,
					normalized_point_src,
					normalized_point_dst;

				x0 = *(points_ptr++);
				y0 = *(points_ptr++);
				x1 = *(points_ptr++);
				y1 = *(points_ptr++);

				point_src << x0, y0, 1.0; // Homogeneous point in the first image
				point_dst << x1, y1, 1.0; // Homogeneous point in the second image

				// Normalized homogeneous point in the first image
				normalized_point_src =
					inverse_intrinsics_src * point_src;
				// Normalized homogeneous point in the second image
				normalized_point_dst =
					inverse_intrinsics_dst * point_dst;

				// The second four columns contain the normalized coordinates.
				*(normalized_points_ptr++) = normalized_point_src(0);
				*(normalized_points_ptr++) = normalized_point_src(1);
				*(normalized_points_ptr++) = normalized_point_dst(0);
				*(normalized_points_ptr++) = normalized_point_dst(1);
			}
		}

		void normalizeImagePoints(const cv::Mat &points_,
			const Eigen::Matrix3d &intrinsics_,
			cv::Mat &normalized_points_)
		{
			const Eigen::Matrix3d inverse_intrinsics = intrinsics_.inverse();

			// Most likely, this is not the fastest solution, but it does
			// not affect the speed of Graph-cut RANSAC, so not a crucial part of
			// this example.
			double x0, y0, x1, y1;
			for (auto r = 0; r < points_.rows; ++r)
			{
				Eigen::Vector3d point,
					normalized_point;

				x0 = points_.at<double>(r, 0);
				y0 = points_.at<double>(r, 1);

				point << x0, y0, 1.0; // Homogeneous point in the image

				// Normalized homogeneous point in the first image
				normalized_point =
					inverse_intrinsics * point;

				// The second four columns contain the normalized coordinates.
				normalized_points_.at<double>(r, 0) = normalized_point(0);
				normalized_points_.at<double>(r, 1) = normalized_point(1);
			}
		}

		void saveStatisticsToFile(
			const gcransac::utils::RANSACStatistics &statistics_,
			const std::string &source_path_,
			const std::string &destination_path_,
			const std::string &filename_,
#ifdef _WIN32
			const int mode_
#else
			const std::ios_base::openmode mode_
#endif
		)
		{
			std::ofstream file(filename_, mode_);

			if (!file.is_open())
			{
				fprintf(stderr, "A problem occured when saving the statistics to file '%s'.\n", filename_.c_str());
				return;
			}

			file << source_path_ << ";"
				<< destination_path_ << ";"
				<< statistics_.iteration_number << ";"
				<< statistics_.processing_time << ";"
				<< statistics_.inliers.size() << ";"
				<< statistics_.main_sampler_name << ";"
				<< statistics_.local_optimizer_sampler_name
				<< std::endl;
			file.close();
		}

		// Returns the sign of a number
		template <typename T>
		inline int signOfNumber(const T val) {
			return (T(0) < val) - (val < T(0));
		}

		// Calculate the rotation matrix from a given homography, translation,
		// plane normal, and plane distance.
		inline Eigen::Matrix3d computeHomographyRotation(
			const Eigen::Matrix3d& H_normalized,
			const Eigen::Vector3d& tstar,
			const Eigen::Vector3d& n,
			const double v) {
			return H_normalized * (Eigen::Matrix3d::Identity() - (2.0 / v) * tstar * n.transpose());
		}

		inline double computeOppositeOfMinor(
			const Eigen::Matrix3d& matrix_,
			const size_t row_,
			const size_t col_) {
			const size_t col1 = col_ == 0 ? 1 : 0;
			const size_t col2 = col_ == 2 ? 1 : 2;
			const size_t row1 = row_ == 0 ? 1 : 0;
			const size_t row2 = row_ == 2 ? 1 : 2;
			return (matrix_(row1, col2) * matrix_(row2, col1) - matrix_(row1, col1) * matrix_(row2, col2));
		}

		double rotationError(const Eigen::Matrix3d &reference_rotation_,
			const Eigen::Matrix3d &estimated_rotation_)
		{
			constexpr double radian_to_degree_multiplier = 180.0 / M_PI;

			const double trace_R_est_times_R_ref =
				(estimated_rotation_ * reference_rotation_.transpose()).trace();

			double error_cos = 0.5 * (trace_R_est_times_R_ref - 1.0);

			// Avoid invalid values due to numerical errors.
			error_cos = std::clamp(error_cos, -1.0, 1.0);

			return radian_to_degree_multiplier * std::acos(error_cos);
		}

		// Decomposing a given homography matrix into the possible pose parameters
		inline void decomposeHomographyMatrix(
			const Eigen::Matrix3d& homography_, // The homography matrix to be decomposed
			const Eigen::Matrix3d& intrinsics_src_, // The intrinsic parameters of the first camera
			const Eigen::Matrix3d& intrinsics_dst_, // The intrinsic parameters of the second camera
			std::vector<Eigen::Matrix3d>& Rs_dst_src_, // The possible rotation matrices
			std::vector<Eigen::Vector3d>& ts_dst_src_, // The possible translation vectors
			std::vector<Eigen::Vector3d>& normals_) // The possible plane normals
		{
			// Remove calibration from homography.
			Eigen::Matrix3d H_normalized = intrinsics_dst_.inverse() * homography_ * intrinsics_src_;

			// Remove scale from normalized homography.
			Eigen::JacobiSVD<Eigen::Matrix3d> hmatrix_norm_svd(H_normalized);
			H_normalized.array() /= hmatrix_norm_svd.singularValues()[1];

			const Eigen::Matrix3d S = H_normalized.transpose() * H_normalized - Eigen::Matrix3d::Identity();

			// Check if H is rotation matrix.
			const double kMinInfinityNorm = 1e-3;
			if (S.lpNorm<Eigen::Infinity>() < kMinInfinityNorm) {
				Rs_dst_src_.emplace_back(H_normalized);
				ts_dst_src_.emplace_back(Eigen::Vector3d::Zero());
				normals_.emplace_back(Eigen::Vector3d::Zero());
				return;
			}

			const double M00 = computeOppositeOfMinor(S, 0, 0);
			const double M11 = computeOppositeOfMinor(S, 1, 1);
			const double M22 = computeOppositeOfMinor(S, 2, 2);

			const double rtM00 = std::sqrt(M00);
			const double rtM11 = std::sqrt(M11);
			const double rtM22 = std::sqrt(M22);

			const double M01 = computeOppositeOfMinor(S, 0, 1);
			const double M12 = computeOppositeOfMinor(S, 1, 2);
			const double M02 = computeOppositeOfMinor(S, 0, 2);

			const int e12 = signOfNumber(M12);
			const int e02 = signOfNumber(M02);
			const int e01 = signOfNumber(M01);

			const double nS00 = std::abs(S(0, 0));
			const double nS11 = std::abs(S(1, 1));
			const double nS22 = std::abs(S(2, 2));

			const std::array<double, 3> nS{ {nS00, nS11, nS22} };
			const size_t idx = std::distance(nS.begin(), std::max_element(nS.begin(), nS.end()));

			Eigen::Vector3d np1;
			Eigen::Vector3d np2;
			if (idx == 0) {
				np1[0] = S(0, 0);
				np2[0] = S(0, 0);
				np1[1] = S(0, 1) + rtM22;
				np2[1] = S(0, 1) - rtM22;
				np1[2] = S(0, 2) + e12 * rtM11;
				np2[2] = S(0, 2) - e12 * rtM11;
			}
			else if (idx == 1) {
				np1[0] = S(0, 1) + rtM22;
				np2[0] = S(0, 1) - rtM22;
				np1[1] = S(1, 1);
				np2[1] = S(1, 1);
				np1[2] = S(1, 2) - e02 * rtM00;
				np2[2] = S(1, 2) + e02 * rtM00;
			}
			else if (idx == 2) {
				np1[0] = S(0, 2) + e01 * rtM11;
				np2[0] = S(0, 2) - e01 * rtM11;
				np1[1] = S(1, 2) + rtM00;
				np2[1] = S(1, 2) - rtM00;
				np1[2] = S(2, 2);
				np2[2] = S(2, 2);
			}

			const double traceS = S.trace();
			const double v = 2.0 * std::sqrt(1.0 + traceS - M00 - M11 - M22);

			const double ESii = signOfNumber(S(idx, idx));
			const double r_2 = 2 + traceS + v;
			const double nt_2 = 2 + traceS - v;

			const double r = std::sqrt(r_2);
			const double n_t = std::sqrt(nt_2);

			const Eigen::Vector3d n1 = np1.normalized();
			const Eigen::Vector3d n2 = np2.normalized();

			const double half_nt = 0.5 * n_t;
			const double esii_t_r = ESii * r;

			const Eigen::Vector3d t1_star = half_nt * (esii_t_r * n2 - n_t * n1);
			const Eigen::Vector3d t2_star = half_nt * (esii_t_r * n1 - n_t * n2);

			const Eigen::Matrix3d R1 = computeHomographyRotation(H_normalized, t1_star, n1, v);
			const Eigen::Vector3d t1 = R1 * t1_star;

			const Eigen::Matrix3d R2 = computeHomographyRotation(H_normalized, t2_star, n2, v);
			const Eigen::Vector3d t2 = R2 * t2_star;

			Rs_dst_src_.emplace_back(-R1);
			Rs_dst_src_.emplace_back(R1);
			Rs_dst_src_.emplace_back(-R2);
			Rs_dst_src_.emplace_back(R2);

			ts_dst_src_.emplace_back(t1);
			ts_dst_src_.emplace_back(-t1);
			ts_dst_src_.emplace_back(t2);
			ts_dst_src_.emplace_back(-t2);

			normals_.emplace_back(-n1);
			normals_.emplace_back(n1);
			normals_.emplace_back(-n2);
			normals_.emplace_back(n2);
		}
		
		inline void decomposeEssentialMatrix(
			const Eigen::Matrix3d& essential_matrix_,
			Eigen::Matrix3d& rotation_1_,
			Eigen::Matrix3d& rotation_2_,
			Eigen::Vector3d& translation_)
		{
			Eigen::Matrix3d d;
			d << 0, 1, 0, -1, 0, 0, 0, 0, 1;

			const Eigen::JacobiSVD<Eigen::Matrix3d> svd(
				essential_matrix_, Eigen::ComputeFullU | Eigen::ComputeFullV);
			Eigen::Matrix3d U = svd.matrixU();
			Eigen::Matrix3d V = svd.matrixV();
			if (U.determinant() < 0) {
				U.col(2) *= -1.0;
			}

			if (V.determinant() < 0) {
				V.col(2) *= -1.0;
			}

			// Possible configurations.
			rotation_1_ = U * d * V.transpose();
			rotation_2_ = U * d.transpose() * V.transpose();
			translation_ = U.col(2).normalized();
		}
	}
}

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    
    if (nrhs != 1) {
        mexErrMsgTxt("One input argument required.");
    }
    if (nlhs > 3) {
        mexErrMsgTxt("Too many output arguments.");
    }
    double const * const matrixAux= static_cast<double const *>(mxGetData(prhs[0]));
    const mwSize *sizeInputMatrix= mxGetDimensions(prhs[0]);
    double*  xX= (double*)malloc(sizeInputMatrix[0] *sizeInputMatrix[1]* sizeof(double));
    double *pindR, *pindt;
    const int size0 = sizeInputMatrix[0]; // Const makes compiler optimization kick in
    const int size1 = sizeInputMatrix[1];
    
    for (int k = 0; k < size0; k++)
    {
        int kOffset = k*size1; // this saves re-computation time
        for (int i = 0; i < size1; i++)
        {
            int iOffset = i*size0;
            xX[i  + kOffset] = matrixAux[i  + kOffset];
        }
    }
    pose_lib::CameraPoseVector solutions;
    
    Eigen::Matrix3d H, K1, K2;
    H << xX[0],xX[1],xX[2],xX[3],xX[4],xX[5],xX[6],xX[7],xX[8];
    K1 << 1,0,0,0,1,0,0,0,1;
    K2 << 1,0,0,0,1,0,0,0,1;
    
    std::vector<Eigen::Matrix3d> Rs_dst_src_;
    x1.push_back(x11); x1.push_back(x12); x1.push_back(x13);x1.push_back(x14); x1.push_back(x15);
    std::vector<Eigen::Vector3d> ts_dst_src_, normals_;
    x2.push_back(x21); x2.push_back(x22); x2.push_back(x23);x2.push_back(x24); x2.push_back(x25);
    
    int noOfSol ;
    
    std::vector<Eigen::Matrix3d> essential_matrices;
    
    decomposeHomographyMatrix(H, // The homography matrix to be decomposed
			K1, // The intrinsic parameters of the first camera
			K2, // The intrinsic parameters of the second camera
			Rs_dst_src_, // The possible rotation matrices
			ts_dst_src_, // The possible translation vectors
			normals_) // The possible plane normals
            
    plhs[0] = mxCreateDoubleMatrix(3, noOfSol*3, mxREAL);
    plhs[1] = mxCreateDoubleMatrix(3, noOfSol, mxREAL);
    
    pindR = mxGetPr(plhs[0]);
    pindt = mxGetPr(plhs[1]);
    
    
//     for(int i =0; i<noOfSol; i++) {
//         for (int k = 0; k < 3; k++) {
//             for (int j = 0; j < 3; j++) {
//                 pindR[(k)+3*j+9*i] = essential_matrices[i](k,j);
//             }
//         }
//     }
    
}