#include <cmath>
#include <iostream>
#include <fstream>
#include <vector>
#include <ctime>
#include <algorithm>

#include <pcl/common/transforms.h>   

#include <pcl/io/pcd_io.h>
#include <pcl/io/vtk_lib_io.h>

#include <pcl/point_cloud.h>

#include <pcl/search/kdtree.h>

#include <pcl/octree/octree.h>

#include <pcl/features/normal_3d_omp.h>
#include <pcl/features/boundary.h>

#include <boost/thread/thread.hpp>

#include <pcl/visualization/pcl_visualizer.h>

#include <pcl/segmentation/region_growing.h>
#include <pcl/segmentation/sac_segmentation.h>
#include <pcl/segmentation/extract_clusters.h>

#include <pcl/filters/statistical_outlier_removal.h>
#include <pcl/filters/radius_outlier_removal.h>
#include <pcl/filters/passthrough.h>
#include <pcl/filters/voxel_grid.h>

#include <pcl/surface/on_nurbs/fitting_curve_2d.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_pdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_tdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_sdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_apdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_asdm.h>
#include <pcl/surface/on_nurbs/fitting_curve_2d_atdm.h>
#include <pcl/surface/on_nurbs/triangulation.h>

using namespace std;

void vectorToUReuler(pcl::PointCloud<pcl::PointXYZINormal>::Ptr cloudNormal)
{
	Eigen::Vector3d dirvec;
	Eigen::Vector3d rpy;
	double roll, pitch = 0, yaw = 0;
	for (size_t i = 0; i < cloudNormal->size(); i++)
	{
		dirvec(0) = cloudNormal->points[i].normal_x;
		dirvec(1) = cloudNormal->points[i].normal_y;
		dirvec(2) = cloudNormal->points[i].normal_z;
		dirvec = dirvec.normalized();
		roll = acos(dirvec.z());
		if (dirvec.x() == 0 && dirvec.y() == 0)
		{
			yaw = 0;
		}
		else
		{
			Eigen::Vector3d noZ = Eigen::Vector3d{ dirvec.x(), dirvec.y(), 0.0 };
			noZ = noZ.normalized();
			if (dirvec.x() <= 0)
				yaw = -acos(-noZ.y());
			else
				yaw = acos(-noZ.y());
		}
		if (pitch < -M_PI)
			pitch = 2 * M_PI + pitch;
		else if (pitch > M_PI)
			pitch = pitch - 2 * M_PI;
		////////////////////////////
		rpy = { roll, pitch, yaw };
		if (roll == 0 && pitch == 0 && yaw == 0)
			rpy = { 0.0, 0.0, 0.0 };
		Eigen::Matrix3d RollM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d YawM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d PitchM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d PRM = Eigen::Matrix3d::Identity();
		RollM(0, 0) = 1.0;
		RollM(0, 1) = 0.0;
		RollM(0, 2) = 0.0;
		RollM(1, 0) = 0.0;
		RollM(1, 1) = cos(roll);
		RollM(1, 2) = -sin(roll);
		RollM(2, 0) = 0.0;
		RollM(2, 1) = sin(roll);
		RollM(2, 2) = cos(roll);

		PitchM(0, 0) = cos(pitch);
		PitchM(0, 1) = 0.0;
		PitchM(0, 2) = sin(pitch);
		PitchM(1, 0) = 0.0;
		PitchM(1, 1) = 1.0;
		PitchM(1, 2) = 0.0;
		PitchM(2, 0) = -sin(pitch);
		PitchM(2, 1) = 0.0;
		PitchM(2, 2) = cos(pitch);

		YawM(0, 0) = cos(yaw);
		YawM(0, 1) = -sin(yaw);;
		YawM(0, 2) = 0.0;
		YawM(1, 0) = sin(yaw);
		YawM(1, 1) = cos(yaw);
		YawM(1, 2) = 0.0;
		YawM(2, 0) = 0.0;
		YawM(2, 1) = 0.0;
		YawM(2, 2) = 1.0;

		Rot = YawM * PitchM * RollM;

		double rotSum = Rot(0, 0) + Rot(1, 1) + Rot(2, 2) - 1;
		double alpha = acos(rotSum / 2);
		double theta = 0.0;
		if (roll >= 0)
			theta = alpha;
		else
			theta = 2.0 * M_PI - alpha;
		double my = 1.0 / (2.0 * sin(theta));

		double rx = my * (Rot(2, 1) - Rot(1, 2)) * theta;
		double ry = my * (Rot(0, 2) - Rot(2, 0)) * theta;
		double rz = my * (Rot(1, 0) - Rot(0, 1)) * theta;
		cout << "UR rad :" << rx << " " << ry << " " << rz << " RPY :" << rpy.x() << " " << rpy.y() << " " << rpy.z() << endl;
	}
}
void vectorToUReuler(pcl::PointXYZ point)
{
	Eigen::Vector3d dirvec;
	Eigen::Vector3d rpy;
	double roll, pitch = 0, yaw = 0;
	{
		dirvec(0) = point.x;
		dirvec(1) = point.y;
		dirvec(2) = point.z;
		dirvec = dirvec.normalized();
		roll = acos(dirvec.z());
		if (dirvec.x() == 0 && dirvec.y() == 0)
		{
			yaw = 0;
		}
		else
		{
			Eigen::Vector3d noZ = Eigen::Vector3d{ dirvec.x(), dirvec.y(), 0.0 };
			noZ = noZ.normalized();
			if (dirvec.x() <= 0)
				yaw = -acos(-noZ.y());
			else
				yaw = acos(-noZ.y());
		}
		if (pitch < -M_PI)
			pitch = 2 * M_PI + pitch;
		else if (pitch > M_PI)
			pitch = pitch - 2 * M_PI;
		////////////////////////////
		rpy = { roll, pitch, yaw };
		if (roll == 0 && pitch == 0 && yaw == 0)
			rpy = { 0.0, 0.0, 0.0 };
		Eigen::Matrix3d RollM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d YawM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d PitchM = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d Rot = Eigen::Matrix3d::Identity();
		Eigen::Matrix3d PRM = Eigen::Matrix3d::Identity();
		RollM(0, 0) = 1.0;
		RollM(0, 1) = 0.0;
		RollM(0, 2) = 0.0;
		RollM(1, 0) = 0.0;
		RollM(1, 1) = cos(roll);
		RollM(1, 2) = -sin(roll);
		RollM(2, 0) = 0.0;
		RollM(2, 1) = sin(roll);
		RollM(2, 2) = cos(roll);

		PitchM(0, 0) = cos(pitch);
		PitchM(0, 1) = 0.0;
		PitchM(0, 2) = sin(pitch);
		PitchM(1, 0) = 0.0;
		PitchM(1, 1) = 1.0;
		PitchM(1, 2) = 0.0;
		PitchM(2, 0) = -sin(pitch);
		PitchM(2, 1) = 0.0;
		PitchM(2, 2) = cos(pitch);

		YawM(0, 0) = cos(yaw);
		YawM(0, 1) = -sin(yaw);;
		YawM(0, 2) = 0.0;
		YawM(1, 0) = sin(yaw);
		YawM(1, 1) = cos(yaw);
		YawM(1, 2) = 0.0;
		YawM(2, 0) = 0.0;
		YawM(2, 1) = 0.0;
		YawM(2, 2) = 1.0;

		Rot = YawM * PitchM * RollM;

		double rotSum = Rot(0, 0) + Rot(1, 1) + Rot(2, 2) - 1;
		double alpha = acos(rotSum / 2);
		double theta = 0.0;
		if (roll >= 0)
			theta = alpha;
		else
			theta = 2.0 * M_PI - alpha;
		double my = 1.0 / (2.0 * sin(theta));

		double rx = my * (Rot(2, 1) - Rot(1, 2)) * theta;
		double ry = my * (Rot(0, 2) - Rot(2, 0)) * theta;
		double rz = my * (Rot(1, 0) - Rot(0, 1)) * theta;
		cout << "UR rad :" << rx << " " << ry << " " << rz << endl;
		//cout << " RPY :" << rpy.x() << " " << rpy.y() << " " << rpy.z() << endl;
	}
}
void STLToPCD()
{
	vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	reader->SetFileName("data...\\data.stl");
	reader->Update();
	vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	polydata = reader->GetOutput();
	polydata->GetNumberOfPoints();

	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>());
	//从ply转pcd
	pcl::io::vtkPolyDataToPointCloud(polydata, *cloud);
	pcl::io::savePCDFileASCII("data...\\data.pcd", *cloud);
}

pcl::PointCloud<pcl::PointXYZ>::Ptr Downsampling(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, float LeafSize = 5.0)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_filtered(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::VoxelGrid<pcl::PointXYZ> vg;
	vg.setInputCloud(SourceCloud);
	vg.setLeafSize(LeafSize, LeafSize, LeafSize);
	vg.filter(*cloud_filtered);
	std::cout << "PointCloud after filtering has: " << cloud_filtered->points.size() << " data points." << std::endl;
	return cloud_filtered;
}

pcl::PointCloud<pcl::PointXYZINormal>::Ptr FindNormal(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, int Radius = 7, int thread = 2)
{
	pcl::PointCloud<pcl::PointXYZINormal>::Ptr cloudNormal(new pcl::PointCloud<pcl::PointXYZINormal>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr filtered(new pcl::PointCloud< pcl::PointXYZ>);
	pcl::PointCloud<pcl::Normal>::Ptr Normals(new pcl::PointCloud<pcl::Normal>);
	pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> nmlOMP;
	nmlOMP.setNumberOfThreads(thread);
	nmlOMP.setInputCloud(SourceCloud);
	nmlOMP.setRadiusSearch(Radius);
	//Eigen::Vector4f centroid;
	//pcl::compute3DCentroid(*SourceCloud, centroid);
	//nmlOMP.setViewPoint(centroid[0], centroid[1], centroid[2]);
	nmlOMP.compute(*Normals);
	pcl::concatenateFields(*SourceCloud, *Normals, *cloudNormal);
	return cloudNormal;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr EdgeDetection(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, int Accurucy = 50.0)
{
	std::cout << "points sieze is:" << SourceCloud->size() << std::endl;
	pcl::PointCloud<pcl::Boundary> boundaries;
	pcl::PointCloud<pcl::Normal>::Ptr normals(new pcl::PointCloud<pcl::Normal>);
	pcl::BoundaryEstimation<pcl::PointXYZ, pcl::Normal, pcl::Boundary> est;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>());
	/*
	pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;  //创建一个快速k近邻查询,查询的时候若该点在点云中，则第一个近邻点是其本身
	kdtree.setInputCloud(cloud);
	int k =2;
	float everagedistance =0;
	for (int i =0; i < cloud->size()/2;i++)
	{
			vector<int> nnh ;
			vector<float> squaredistance;
			//  pcl::PointXYZ p;
			//   p = cloud->points[i];
			kdtree.nearestKSearch(cloud->points[i],k,nnh,squaredistance);
			everagedistance += sqrt(squaredistance[1]);
			//   cout<<everagedistance<<endl;
	}
	everagedistance = everagedistance/(cloud->size()/2);
	cout<<"everage distance is : "<<everagedistance<<endl;
*/
	pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normEst;  //其中pcl::PointXYZ表示输入类型数据，pcl::Normal表示输出类型,且pcl::Normal前三项是法向，最后一项是曲率
	normEst.setInputCloud(SourceCloud);
	normEst.setSearchMethod(tree);
	// normEst.setRadiusSearch(2);  //法向估计的半径
	normEst.setKSearch(9);  //法向估计的点数
	normEst.compute(*normals);
	cout << "normal size is " << normals->size() << endl;

	est.setInputCloud(SourceCloud);
	est.setInputNormals(normals);
	//  est.setAngleThreshold(90);
	//   est.setSearchMethod (pcl::search::KdTree<pcl::PointXYZ>::Ptr (new pcl::search::KdTree<pcl::PointXYZ>));
	est.setSearchMethod(tree);
	est.setKSearch(Accurucy);  //一般这里的数值越高，最终边界识别的精度越好
	est.compute(boundaries);

	//  pcl::PointCloud<pcl::PointXYZ> boundPoints;
	pcl::PointCloud<pcl::PointXYZ>::Ptr boundPoints(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ> noBoundPoints;
	int countBoundaries = 0;
	for (int i = 0; i < SourceCloud->size(); i++) {
		uint8_t x = (boundaries.points[i].boundary_point);
		int a = static_cast<int>(x); //该函数的功能是强制类型转换
		if (a == 1)
		{
			//  boundPoints.push_back(cloud->points[i]);
			(*boundPoints).push_back(SourceCloud->points[i]);
			countBoundaries++;
		}
		else
			noBoundPoints.push_back(SourceCloud->points[i]);

	}
	std::cout << "boudary size is：" << countBoundaries << std::endl;
	return boundPoints;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr SegEuclidean(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, double threshold = 10.0, int MinClusterSize = 50, int MaxClusterSize = 1500)
{
	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(SourceCloud);
	ec.setClusterTolerance(threshold); // 设置临近搜索的搜索半径
	ec.setMinClusterSize(MinClusterSize);    // 每个簇（集群）的最小大小
	ec.setMaxClusterSize(MaxClusterSize);  // 每个簇（集群）的最大大小
	ec.setSearchMethod(tree);     // 设置点云搜索算法
	ec.setInputCloud(SourceCloud);   // 设置输入点云
	ec.extract(cluster_indices);
	cout << cluster_indices.size() << endl;
	// 每次创建一个新的点云数据集，并且将所有当前簇的点写入到点云数据集中。
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_cluster(new pcl::PointCloud<pcl::PointXYZ>);

	const std::vector<int>& indices = cluster_indices.begin()->indices;
	if (cluster_indices.size() != 0)
	{
		for (std::vector<int>::const_iterator pit = indices.begin(); pit != indices.end(); ++pit)
			cloud_cluster->points.push_back(SourceCloud->points[*pit]);
	}

	return cloud_cluster;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr SegEuclidean_Max(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, float threshold)
{
	pcl::search::KdTree<pcl::PointXYZ>::Ptr tree(new pcl::search::KdTree<pcl::PointXYZ>);
	tree->setInputCloud(SourceCloud);

	std::vector<pcl::PointIndices> cluster_indices;
	pcl::EuclideanClusterExtraction<pcl::PointXYZ> ec;
	ec.setClusterTolerance(threshold); // 2cm
	ec.setMinClusterSize(100);
	ec.setMaxClusterSize(25000);
	ec.setSearchMethod(tree);
	ec.setInputCloud(SourceCloud);
	ec.extract(cluster_indices);

	pcl::PointCloud<pcl::PointXYZ>::Ptr cluster(new pcl::PointCloud<pcl::PointXYZ>);
	std::vector<int> max_cluster;
	int max_size = 0;
	std::cout << cluster_indices.size() << std::endl;
	if (cluster_indices.size() != 0)
	{
		for (std::vector<pcl::PointIndices>::const_iterator it = cluster_indices.begin(); it != cluster_indices.end(); ++it)
		{
			//std::cout << it->indices.size() << std::endl;
			if (it->indices.size() > max_size)
			{
				//std::cout << it->indices.size() << std::endl;
				max_size = it->indices.size();
				max_cluster = it->indices;
			}
		}
	}

	for (std::vector<int>::const_iterator pit = max_cluster.begin(); pit != max_cluster.end(); ++pit)
		cluster->points.push_back(SourceCloud->points[*pit]);

	return cluster;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr RadiusDenoise(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, float Radious, int MinPoint)
{
	//-------------------去噪-------------------//
	pcl::PointCloud<pcl::PointXYZ>::Ptr denoised(new pcl::PointCloud<pcl::PointXYZ>);

	pcl::RadiusOutlierRemoval<pcl::PointXYZ> pcFilter;  //创建滤波器对象
	pcFilter.setInputCloud(SourceCloud);             //设置待滤波的点云
	pcFilter.setRadiusSearch(Radious);               // 设置搜索半径
	pcFilter.setMinNeighborsInRadius(MinPoint);      // 设置一个内点最少的邻居数目

	pcFilter.filter(*denoised);
	return denoised;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr MeanHeight(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, int MeanPointNum, bool isArrange, bool Xsurface = false, bool Ysurface = false, bool Zsurface = true)
{
	double tmp;
	pcl::PointXYZ point;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);
	if (!isArrange)
	{
		pcl::KdTreeFLANN<pcl::PointXYZ> kdtree;
		kdtree.setInputCloud(SourceCloud);

		int K = MeanPointNum;
		std::vector<int> pointIdxNKNSearch(K);
		std::vector<float> pointNKNSquaredDistance(K);

		for (size_t i = 0; i < SourceCloud->size(); i++)
		{
			tmp = 0;
			kdtree.nearestKSearch(SourceCloud->points[i], K, pointIdxNKNSearch, pointNKNSquaredDistance);
			for (size_t j = 0; j < pointIdxNKNSearch.size(); ++j)
			{
				if (Xsurface)
				{
					tmp += SourceCloud->points[pointIdxNKNSearch[j]].x;
				}
				if (Ysurface)
				{
					tmp += SourceCloud->points[pointIdxNKNSearch[j]].y;
				}
				if (Zsurface)
				{
					tmp += SourceCloud->points[pointIdxNKNSearch[j]].z;
				}
			}
			tmp /= pointIdxNKNSearch.size();
			if (Xsurface)
			{
				point.x = tmp;
				point.y = SourceCloud->points[i].y;
				point.z = SourceCloud->points[i].z;
			}
			if (Ysurface)
			{
				point.x = SourceCloud->points[i].x;
				point.y = tmp;
				point.z = SourceCloud->points[i].z;
			}
			if (Zsurface)
			{
				point.x = SourceCloud->points[i].x;
				point.y = SourceCloud->points[i].y;
				point.z = tmp;
			}
			cloud_After->points.push_back(point);
		}
	}

	else
	{
		int count = SourceCloud->points.size();
		int MinRange, MaxRange, Total;
		double tmpHeight;
		Total = 2 * MeanPointNum + 1;

		for (int i = 0; i < count; i++)
		{
			tmpHeight = 0;
			for (int k = i - MeanPointNum; k < i + MeanPointNum + 1; k++) {

				int index = k;
				if (k < 0) {
					index = count + k;
				}
				else if (k >= count) {
					index = k - count;
				}
				tmpHeight += SourceCloud->points[index].z;
			}
			tmpHeight /= Total;
			point.x = SourceCloud->points[i].x;
			point.y = SourceCloud->points[i].y;
			point.z = tmpHeight;
			cloud_After->points.push_back(point);
		}
	}

	return cloud_After;

}

pcl::PointCloud<pcl::PointXYZ>::Ptr AxisRotate(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, float Degree, string Axis = "X")
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);

	Eigen::Matrix4f transform = Eigen::Matrix4f::Identity();
	float theta = (M_PI / 180) * Degree; // The angle of rotation in radians
	if (Axis == "X")
	{
		transform(1, 1) = cos(theta);
		transform(2, 1) = -sin(theta);
		transform(1, 2) = sin(theta);
		transform(2, 2) = cos(theta);
	}
	if (Axis == "Y")
	{
		transform(0, 0) = cos(theta);
		transform(0, 2) = sin(theta);
		transform(2, 0) = -sin(theta);
		transform(2, 2) = cos(theta);
	}
	if (Axis == "Z")
	{
		transform(0, 0) = cos(theta);
		transform(0, 1) = -sin(theta);
		transform(1, 0) = sin(theta);
		transform(1, 1) = cos(theta);
	}

	pcl::transformPointCloud(*SourceCloud, *cloud_After, transform);
	return cloud_After;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr PolarSortingXY(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Process(new pcl::PointCloud<pcl::PointXYZ>);

	pcl::copyPointCloud(*SourceCloud, *cloud_Process);

	for (size_t i = 0; i < cloud_Process->points.size(); i++)
	{
		cloud_Process->points[i].z = 0;
	}

	Eigen::Vector4d centroid;

	pcl::compute3DCentroid(*cloud_Process, centroid);

	std::cout << "The XYZ coordinates of the centroid are: ("
		<< centroid[0] << ", "
		<< centroid[1] << ", "
		<< centroid[2] << ")." << std::endl;

	double x, y, tmp;
	double* theta;
	theta = new double[cloud_Process->points.size()];

	for (size_t i = 0; i < cloud_Process->points.size(); i++)
	{
		x = (double)cloud_Process->points[i].x - centroid[0];
		y = (double)cloud_Process->points[i].y - centroid[1];
		theta[i] = atan2(y, x);
	}

	sort(theta, theta + cloud_Process->points.size());

	for (size_t i = 0; i < cloud_Process->points.size(); i++)
	{
		for (size_t j = 0; j < SourceCloud->points.size(); j++)
		{
			x = (double)SourceCloud->points[j].x - centroid[0];
			y = (double)SourceCloud->points[j].y - centroid[1];
			tmp = theta[i] - atan2(y, x);
			//cout << theta[i] << " - " << atan2(y, x) << " = " << tmp << endl;
			if (tmp == 0)
			{
				cloud_After->points.push_back(SourceCloud->points[j]);
				break;
			}
		}
	}

	delete[] theta;
	cloud_After->height = 1;
	cloud_After->width = cloud_After->points.size();

	return cloud_After;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr ArrangeRouteXY_Deg(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, double DegreeGap)
{
	pcl::PointXYZ point;
	Eigen::Vector4d centroid;
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Process(new pcl::PointCloud<pcl::PointXYZ>);

	pcl::copyPointCloud(*SourceCloud, *cloud_Process);
	for (size_t i = 0; i < cloud_Process->points.size(); i++)
	{
		cloud_Process->points[i].z = 0;
	}

	pcl::compute3DCentroid(*cloud_Process, centroid);
	point.x = centroid[0];
	point.y = centroid[1];
	point.z = centroid[2];

	double theta, x, y;

	for (size_t deg = 0; deg < 360; deg += DegreeGap)
	{
		for (size_t i = 0; i < cloud_Process->points.size(); i++)
		{
			x = (double)cloud_Process->points[i].x - centroid[0];
			y = (double)cloud_Process->points[i].y - centroid[1];
			theta = atan2(y, x) * 180 / M_PI;

			if (theta < 0)
			{
				theta += 360;
			}

			if (abs(deg - theta) < 0.5)
			{
				cloud_After->points.push_back(SourceCloud->points[i]);
				break;
			}
		}
	}

	return cloud_After;
}

pcl::PointCloud<pcl::PointXYZ>::Ptr ArrangeRouteXY_Ratio(pcl::PointCloud<pcl::PointXYZ>::Ptr SourceCloud, double TotalPoint)
{
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);

	int gap = SourceCloud->points.size() / TotalPoint;
	if (gap <= 1)
	{
		cloud_After = SourceCloud;
	}
	else
	{
		for (size_t i = 0; i < SourceCloud->points.size(); i += gap)
		{
			cloud_After->points.push_back(SourceCloud->points[i]);
			if (cloud_After->points.size() >= TotalPoint + 1)
			{
				break;
			}
		}
	}
	return cloud_After;
}

double rounding(float num, int index)
{
	//https://dotblogs.com.tw/forloop/2016/07/31/rounding
	bool isNegative = false; // whether is negative number or not

	if (num < 0) // if this number is negative, then convert to positive number
	{
		isNegative = true;
		num = -num;
	}

	if (index >= 0)
	{
		int multiplier;
		multiplier = pow(10, index);
		num = (int)(num * multiplier + 0.5) / (multiplier * 1.0);
	}

	if (isNegative) // if this number is negative, then convert to negative number
	{
		num = -num;
	}

	return num;
}

boost::shared_ptr<pcl::visualization::PCLVisualizer> viewer(new pcl::visualization::PCLVisualizer("Cloud Viewer"));
int main(int argc, char** argv)
{
#pragma region Variable
	srand((unsigned int)time(NULL));
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_corner(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_down(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_After(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Edge(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_Centroid(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_down_rotate(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_down_rotate_sort(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr cloud_down_rotate_sort_Y0(new pcl::PointCloud<pcl::PointXYZ>);

	pcl::PointCloud<pcl::PointXYZ>::Ptr edgeDetected(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr denoised(new pcl::PointCloud<pcl::PointXYZ>);
	pcl::PointCloud<pcl::PointXYZ>::Ptr edge(new pcl::PointCloud<pcl::PointXYZ>);

	pcl::PointCloud<pcl::PointXYZ>::Ptr filtered(new pcl::PointCloud< pcl::PointXYZ>);
	pcl::PointCloud<pcl::Normal>::Ptr Normals(new pcl::PointCloud<pcl::Normal>);
	pcl::PointCloud<pcl::PointXYZINormal>::Ptr cloudNormal(new pcl::PointCloud<pcl::PointXYZINormal>);
	pcl::PointCloud<pcl::PointXYZINormal>::Ptr tmpcloudNormal(new pcl::PointCloud<pcl::PointXYZINormal>);
	pcl::PointCloud<pcl::PointXYZINormal>::Ptr tmpcloudNormalA(new pcl::PointCloud<pcl::PointXYZINormal>);

	pcl::PCDReader reader;
#pragma endregion

#pragma region Cube pose rotate
	////-------------------Cube Test-------------------//
	//pcl::PointCloud<pcl::Normal>::Ptr Normals(new pcl::PointCloud<pcl::Normal>);
	//pcl::io::loadPCDFile<pcl::PointXYZ>("D:\\Program\\[3]_PCL\\PCD_Data\\cube200.pcd", *cloud);
	//pcl::NormalEstimationOMP<pcl::PointXYZ, pcl::Normal> nmlOMP;
	//nmlOMP.setNumberOfThreads(2);
	//nmlOMP.setInputCloud(cloud);
	//nmlOMP.setRadiusSearch(7);
	//Eigen::Vector4f centroid;
	//pcl::compute3DCentroid(*filtered, centroid);
	//nmlOMP.setViewPoint(100.0, 100.0, 100.0);
	//nmlOMP.compute(*Normals);
	//pcl::concatenateFields(*cloud, *Normals, *tmpcloudNormal);
	////-------------------CenterPoint-------------------//
	//for (size_t i = 0; i < cloud->size(); i++)
	//{
	//	if ((cloud->points[i].x == 0.0 && cloud->points[i].y == 0.0 && cloud->points[i].z == 0.0) ||
	//		(cloud->points[i].x == 200.0 && cloud->points[i].y == 0.0 && cloud->points[i].z == 0.0) ||
	//		(cloud->points[i].x == 200.0 && cloud->points[i].y == 200.0 && cloud->points[i].z == 0.0) ||
	//		(cloud->points[i].x == 0.0 && cloud->points[i].y == 200.0 && cloud->points[i].z == 0.0))
	//	{
	//		cloudNormal->push_back(tmpcloudNormal->points[i]);
	//	}
	//}
	//pcl::io::savePCDFile("D:\\Program\\[3]_PCL\\PCD_Data\\cube_normal.pcd", *cloudNormal);
	//////-------------------UR Rotate-------------------//

	//pcl::io::loadPCDFile<pcl::PointXYZINormal>("D:\\Program\\[3]_PCL\\PCD_Data\\cube_normal.pcd", *cloudNormal);
	//vectorToUReuler(cloudNormal);
#pragma endregion

#pragma region Calculate corner point
	//reader.read("D:\\JasonWork\\PCD\\Formal\\curve_downsample_leaf3n.pcd", *cloud);
	//reader.read("D:\\JasonWork\\PCD\\Formal\\Edge_Normal.pcd", *cloudNormal);
	//reader.read("D:\\JasonWork\\PCD\\Formal\\EdgeXY.pcd", *cloud_Edge);

	//cloud = AxisRotate(cloud, 90);

	//double xVec, yVec, distY, distX, theta, maxY = 0.0, maxX = 0.0;
	//double degminY = 5.0, degminX = 5.0;
	//int YmaxIndex = 0, XmaxIndex = 0, YminIndex = 0, XminIndex = 0;
	//pcl::PointXYZ minPt, maxPt;
	//pcl::getMinMax3D(*cloud_Edge, minPt, maxPt);

	//for (size_t i = 0; i < cloudNormal->points.size(); i++)
	//{
	//	if (abs(cloudNormal->points[i].x - maxPt.x) < 0.1)
	//	{
	//		tmpcloudNormal->points.push_back(cloudNormal->points[i]);
	//	}
	//	if (abs(cloudNormal->points[i].x - minPt.x) < 0.1)
	//	{
	//		tmpcloudNormal->points.push_back(cloudNormal->points[i]);
	//	}
	//}
	//for (size_t i = 0; i < cloud_Edge->points.size(); i++)
	//{
	//	if (cloud_Edge->points[i].x > 0)
	//	{
	//		xVec = cloud_Edge->points[i].x - tmpcloudNormal->points[0].x;
	//		yVec = cloud_Edge->points[i].y - tmpcloudNormal->points[0].y;
	//		//dist = sqrtf(pow(xVec, 2) + pow(yVec, 2));
	//		theta = atan2(yVec, xVec);
	//		//cout << "Polar Angle : " << theta * 180.0 / M_PI << endl;
	//		if (abs(theta * 180.0 / M_PI) > 114 && abs(theta * 180.0 / M_PI) < 115)
	//		{
	//			distY = abs(cloud_Edge->points[i].y - tmpcloudNormal->points[0].y);
	//			if (distY > maxY)
	//			{
	//				maxY = distY;
	//				YmaxIndex = i;
	//			}
	//			edge->points.push_back(cloud_Edge->points[i]);
	//		}
	//	}
	//	if (cloud_Edge->points[i].x < 0)
	//	{
	//		xVec = cloud_Edge->points[i].x - tmpcloudNormal->points[1].x;
	//		yVec = cloud_Edge->points[i].y - tmpcloudNormal->points[1].y;
	//		//dist = sqrtf(pow(xVec, 2) + pow(yVec, 2));
	//		theta = atan2(yVec, xVec);
	//		//cout << "Polar Angle : " << theta * 180.0 / M_PI << endl;
	//		if (abs(theta * 180.0 / M_PI) > 64 && abs(theta * 180.0 / M_PI) < 65)
	//		{
	//			distX = abs(cloud_Edge->points[i].x - tmpcloudNormal->points[1].x);
	//			if (distX > maxX)
	//			{
	//				maxX = distX;
	//				XmaxIndex = i;
	//			}
	//			edge->points.push_back(cloud_Edge->points[i]);
	//		}
	//	}

	//	if (cloud_Edge->points[i].y > 0)
	//	{

	//		if (abs(cloud_Edge->points[i].x) < degminX)
	//		{
	//			degminX = abs(cloud_Edge->points[i].x);
	//			XminIndex = i;
	//		}
	//	}
	//	if (cloud_Edge->points[i].y < 0)
	//	{

	//		if (abs(cloud_Edge->points[i].x) < degminY)
	//		{
	//			degminY = abs(cloud_Edge->points[i].x);
	//			YminIndex = i;
	//		}
	//	}
	//
	//}
	//tmpcloudNormal->points.push_back(cloudNormal->points[YmaxIndex]);
	//tmpcloudNormal->points.push_back(cloudNormal->points[XmaxIndex]);
	//tmpcloudNormal->points.push_back(cloudNormal->points[YminIndex]);
	//tmpcloudNormal->points.push_back(cloudNormal->points[XminIndex]);
#pragma endregion

#pragma region Create Normal vector and match with Edge points

	reader.read("D:\\JasonWork\\PCD\\Formal\\curve_downsample_leaf3n.pcd", *cloud);
	reader.read("D:\\JasonWork\\PCD\\Formal\\Edge.pcd", *cloud_Edge);
	cloud = AxisRotate(cloud, 90);
	cloud_Edge = AxisRotate(cloud_Edge, 90);
	cloud_Edge = PolarSortingXY(cloud_Edge);
	cloud_Edge = MeanHeight(cloud_Edge, 5, true);

	tmpcloudNormal = FindNormal(cloud, 100);

	cloud_Edge = ArrangeRouteXY_Ratio(cloud_Edge, 20);

	pcl::PointXYZ point;

	double xGate, yGate;

	for (size_t i = 0; i < cloud_Edge->points.size(); i++)
	{
		for (size_t j = 0; j < tmpcloudNormal->points.size(); j++)
		{
			xGate = cloud_Edge->points[i].x - tmpcloudNormal->points[j].x;
			yGate = cloud_Edge->points[i].y - tmpcloudNormal->points[j].y;
			if (abs(xGate) < 0.5 && abs(yGate) < 0.5)
			{
				cloudNormal->points.push_back(tmpcloudNormal->points[j]);
				if (cloudNormal->points[i].normal_z > 0)
				{
					cloudNormal->points[i].normal_z *= -1.0;
				}
				break;
			}
		}
	}

	for (size_t i = 0; i < cloud->points.size(); i++)
	{
		cloud->points[i].x += -140.5;
		cloud->points[i].y += -582.61;
		cloud->points[i].z += 80.7;
	}
	for (size_t i = 0; i < cloudNormal->points.size(); i++)
	{
		cloudNormal->points[i].x += -140.5;
		cloudNormal->points[i].y += -582.61;
		cloudNormal->points[i].z += 80.7;



		cloud_Edge->points[i].x += -140.5;
		cloud_Edge->points[i].y += -582.61;
		cloud_Edge->points[i].z += 80.7;

		cout << "X : " << cloudNormal->points[i].x << endl;
		cout << "Y : " << cloudNormal->points[i].y << endl;
		cout << "Z : " << cloudNormal->points[i].z << endl;
		point.x = cloudNormal->points[i].normal_x;
		point.y = cloudNormal->points[i].normal_y;
		point.z = cloudNormal->points[i].normal_z;
		vectorToUReuler(point);
	}
	//tmpcloudNormal->height = 1;
	//tmpcloudNormal->width = tmpcloudNormal->points.size();
	//pcl::io::savePCDFileASCII("C:\\Users\\JasonYJSu\\Desktop\\PCD\\Formal\\Edge_Normal.pcd", *cloudNormal);
#pragma endregion

#pragma region Transfer2.csv and find the corner point
	//reader.read("C:\\Users\\JasonYJSu\\Desktop\\PCD\\Formal\\Edge.pcd", *cloud);
	//pcl::copyPointCloud(*cloud, *cloud_After);

	//cloud = AxisRotate(cloud, 90);

	//cloud_After = AxisRotate(cloud_After, 90);

	//cloud_After = MeanHeight(cloud_After, 20);

	//pcl::copyPointCloud(*cloud_After, *cloud_down);

	//for (size_t i = 0; i < cloud_down->points.size(); i++)
	//{
	//	cloud_down->points[i].z = 0;
	//}

	//int maxX, maxY;
	//float dist, tmpdist;
	//bool isFind = false;
	//pcl::PointXYZ minPt, maxPt;
	//pcl::getMinMax3D(*cloud_down, minPt, maxPt);
	//for (size_t i = 0; i < cloud_down->points.size(); i++)
	//{
	//	cloud_down->points[i].x -= minPt.x;
	//	cloud_down->points[i].y -= minPt.y;
	//}

	//Eigen::Vector4d centroid;

	//pcl::compute3DCentroid(*cloud_down, centroid);

	//std::cout << "The XYZ coordinates of the centroid are: ("
	//	<< centroid[0] << ", "
	//	<< centroid[1] << ", "
	//	<< centroid[2] << ")." << std::endl;

	//double x, y, tmp;
	//double* theta;
	//double* distance, * dX, * dY, * tmpdistance;
	//theta = new double[cloud_down->points.size()];
	//distance = new double[cloud_down->points.size()];
	//tmpdistance = new double[cloud_down->points.size()];
	//dX = new double[cloud_down->points.size()];
	//dY = new double[cloud_down->points.size()];

	//for (size_t i = 0; i < cloud_down->points.size(); i++)
	//{
	//	x = (double)cloud_down->points[i].x - centroid[0];
	//	y = (double)cloud_down->points[i].y - centroid[1];
	//	theta[i] = atan2(y, x);
	//}

	//sort(theta, theta + cloud_down->points.size());

	//for (size_t i = 0; i < cloud_down->points.size(); i++)
	//{
	//	for (size_t j = 0; j < cloud_down->points.size(); j++)
	//	{
	//		x = (double)cloud_down->points[j].x - centroid[0];
	//		y = (double)cloud_down->points[j].y - centroid[1];
	//		tmp = theta[i] - atan2(y, x);
	//		//cout << theta[i] << " - " << atan2(y, x) << " = " << tmp << endl;
	//		if (tmp == 0)
	//		{
	//			cloud_tmp->points.push_back(cloud_After->points[j]);
	//			break;
	//		}
	//	}
	//}

	//for (size_t i = 0; i < cloud_tmp->points.size(); i++)
	//{
	//	x = centroid[0] - cloud_tmp->points[i].x;
	//	y = centroid[1] - cloud_tmp->points[i].y;
	//	dist = sqrtf(pow(x, 2) + pow(y, 2));

	//	if (i == 112 || i == 835 || i == 1101 || i == 1918)
	//	{
	//		cloud_corner->points.push_back(cloud_tmp->points[i]);
	//	}
	//	dX[i] = cloud_tmp->points[i].x;
	//	dY[i] = cloud_tmp->points[i].y;
	//	tmpdistance[i] = dist;
	//}

	//int count = cloud_tmp->points.size();
	//int OutRange, MaxRange, Total;

	//Total = 2 * POINTLITMIT + 1;

	//for (int i = 0; i < count; i++)
	//{
	//	tmpdist = 0;
	//	for (int k = i - 5; k < i + 6; k++) {

	//		int index = k;
	//		if (k < 0) {
	//			index = count + k;
	//		}
	//		else if (k >= count) {
	//			index = k - count;
	//		}
	//		tmpdist += tmpdistance[index];
	//	}
	//	tmpdist /= Total;
	//	distance[i] = tmpdist;
	//}

	//ofstream ofs;
	//ofs.open("C:\\Users\\JasonYJSu\\Desktop\\sortDistance.txt");
	//for (size_t i = 0; i < count; i++)
	//{
	//	ofs << dX[i] << "," << dY[i] << "," << distance[i] << "\n";
	//}
	//ofs.close();
	//delete[] theta;
	//delete[] distance;
	//delete[] dX;
	//delete[] dY;

#pragma endregion

#pragma region After downsampling original PCD do follow step to find edge
	//////----------------Find Edge----------------//
	////reader.read("C:\\Users\\JasonYJSu\\Desktop\\PCD\\Formal\\curve_downsample_leaf3n.pcd", *cloud_Y0);
	////for (size_t i = 0; i < cloud->size(); i++)
	////{
	////	cloud_Y0->points.push_back(cloud->points[i]);
	////	cloud_Y0->points[i].y = 0;
	////}
	////
	////edgeDetected = EdgeDetection(cloud_Y0);
	////
	////edgeDetected = EdgeDetection(edgeDetected);
	////
	////denoised = RadiusDenoise(edgeDetected, 20.0, 1);
	////
	////edge = SegEuclidean_Max(denoised, 20.0);
	////
	//////cout << cloud->points.size() << endl;
	////
	////float threshValue = 0.5;
	////float xGate, zGate;
	////for (size_t i = 0; i < edge->size(); i++)
	////{
	////	for (size_t j = 0; j < cloud->size(); j++)
	////	{
	////		xGate = abs(edge->points[i].x - cloud->points[j].x);
	////		zGate = abs(edge->points[i].z - cloud->points[j].z);
	////
	////		if (xGate < threshValue && zGate < threshValue)
	////		{
	////			edge->points[i].y = cloud->points[j].y;
	////			//cloud_outedge->points.push_back(cloud_ori->points[j]);
	////			break;
	////		}
	////	}
	////}
	////edge->height = 1;
	////edge->width = edge->points.size();
#pragma endregion

	//////-------------------Save-------------------//

	//pcl::io::savePCDFileASCII("C:\\Users\\JasonYJSu\\Desktop\\PCD\\Formal\\Edge_Sort_40.pcd", *cloud_down_rotate_sort);

	//-------------------Run Time Calculate-------------------//
	//clock_t start, end; // 儲存時間用的變數
	//start = clock(); // 計算開始時間
	//// 這邊放主要計算
	/////////////////////
	//
	/////////////////////
	//end = clock(); // 計算結束時間
	//double cpu_time_used = ((double)(end - start)) / CLOCKS_PER_SEC; // 計算實際花費時間
	//
	//cout << "程式執行時間:" << cpu_time_used << endl;
#pragma region STL to PCD
	////-------------------STL to PCD-------------------//
	//vtkSmartPointer<vtkSTLReader> reader = vtkSmartPointer<vtkSTLReader>::New();
	//reader->SetFileName("C:\\Users\\JasonYJSu\\Desktop\\CAD\\Radiator grille.stl");
	//reader->Update();
	//vtkSmartPointer<vtkPolyData> polydata = vtkSmartPointer<vtkPolyData>::New();
	//polydata = reader->GetOutput();
	//polydata->GetNumberOfPoints();

	//pcl::io::vtkPolyDataToPointCloud(polydata, *cloud);
	//pcl::io::savePCDFileASCII("C:\\Users\\JasonYJSu\\Desktop\\CAD\\curve.pcd", *cloud);
#pragma endregion

#pragma region Segmentation
	////-------------------切割(區域成長)-------------------//
	//pcl::search::Search<pcl::PointXYZ>::Ptr treeseg(new pcl::search::KdTree<pcl::PointXYZ>);
	////求法线　和　曲率　
	//
	//pcl::NormalEstimation<pcl::PointXYZ, pcl::Normal> normal_estimator;
	//normal_estimator.setSearchMethod(treeseg);
	//normal_estimator.setInputCloud(cloud);
	//normal_estimator.setKSearch(7);
	//normal_estimator.compute(*normals);
	//
	////直通滤波在Z轴的0到1米之间 剔除　nan　和　噪点
	//pcl::IndicesPtr indices(new std::vector <int>);//指向int类型的vector类的空智能指针
	//pcl::PassThrough<pcl::PointXYZ> pass;
	//pass.setInputCloud(cloud);
	//pass.setFilterFieldName("x");
	//pass.setFilterLimits(0.0, 1.0);
	//pass.filter(*indices);
	//
	////区域增长聚类分割对象　<点，法线>
	//pcl::RegionGrowing<pcl::PointXYZ, pcl::Normal> reg;
	//reg.setMinClusterSize(100);     //最小的聚类的点数
	//reg.setMaxClusterSize(1000);//最大的聚类的点数
	//reg.setSearchMethod(treeseg);     //搜索方式
	//reg.setNumberOfNeighbours(7); //设置搜索的邻域点的个数
	//reg.setInputCloud(cloud);      //输入点
	//
	//reg.setInputNormals(normals);  //输入的法线
	//reg.setSmoothnessThreshold(3.0 / 180.0 * M_PI);//设置平滑度 法线差值阈值
	////reg.setCurvatureThreshold(1.0);                //设置曲率的阀值
	//
	//std::vector <pcl::PointIndices> clusters;
	//reg.extract(clusters);//提取点的索引
	//
	//std::cout << std::endl;

#pragma endregion

#pragma region Viewer
	////-------------------Viewer-------------------//

	//viewer->setBackgroundColor(0, 0, 0);

	viewer->addPointCloudNormals<pcl::PointXYZINormal>(cloudNormal, 1, -50, "cloudNormal");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloudNormal");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloudNormal");

	//viewer->addPointCloudNormals<pcl::PointXYZINormal>(tmpcloudNormal, 1, -50, "cloudNormal");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "cloudNormal");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "cloudNormal");

	//viewer->addPointCloudNormals<pcl::PointXYZINormal>(tmpcloudNormalA, 1, -50, "tmpcloudNormalA");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.0, 1.0, 0.0, "tmpcloudNormalA");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 5, "tmpcloudNormalA");

	//viewer->addPointCloud(cloud_Centroid, "cloud_Centroid");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloud_Centroid");;
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "cloud_Centroid");

	viewer->addPointCloud(cloud, "cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 0.5, 0.5, 0.5, "cloud");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud");

	//viewer->addPointCloud(cloud_After, "cloud_After");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 1.0, 1.0, "cloud_After");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 2, "cloud_After");

	//viewer->addPointCloud(cloud_corner, "cloud_corner");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "cloud_corner");
	//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 9, "cloud_corner");

	viewer->addPointCloud(cloud_Edge, "edge");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "edge");
	viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_POINT_SIZE, 10, "edge");

	////pcl::PointCloud <pcl::PointXYZRGB>::Ptr colored_cloud = reg.getColoredCloud();
	////viewer->addPointCloud(colored_cloud, "colored_cloud");

	viewer->addCoordinateSystem(100.0);
	viewer->initCameraParameters();
	viewer->setCameraPosition(0, 300, 300, 0, 0, 300);
#pragma endregion
	int a = 0;
	while (!viewer->wasStopped())
	{
		viewer->spinOnce(50);
		//edge->points.push_back(cloud_After->points[a]);
		//viewer->updatePointCloud(edge, "edge");
		//viewer->setPointCloudRenderingProperties(pcl::visualization::PCL_VISUALIZER_COLOR, 1.0, 0.0, 0.0, "edge");
		//a++;
		//edge->points.clear();
		//if (a == cloud_After->points.size())
		//{
		//	a = 0;
		//}
		//boost::this_thread::sleep(boost::posix_time::microseconds(1000));
	}

	return (0);
}
