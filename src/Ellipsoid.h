

#ifndef __ELLIPSOID_INC__
#define __ELLIPSOID_INC__

#include <iostream>
#include <math.h>
#include <Eigen/Dense>


namespace coti {
	class Ellipsoid {
	public:
		Ellipsoid(double a = 6378137.0, double invf = 298.257223563);
		Ellipsoid(const Ellipsoid& ellipsoid);
		Ellipsoid& operator=(const Ellipsoid& ellipsoid);

		double Z_XbarPsi(double xBar, double psi);

		bool GeocentricToGeodetic(double X, double Y, double Z, double& lat, double& lon, double& hei);
		void GeocentricToGeodeticOzone(double X, double Y, double Z, double& lat, double& lon, double& hei);
		void GeodeticToGeocentric(double lat, double lon, double hei, double& X, double& Y, double& Z);

		void ECEFToENU(double X, double Y, double Z, double lat0, double lon0, double h0, double& x, double& y, double& z);
		void ENUToECEF(double x, double y, double z, double lat0, double lon0, double h0, double& X, double& Y, double& Z);

		double a() { return a_; };
		double invf() { return invf_; };

	protected:
		double a_;
		double invf_;
		double f_;
		double b_;
		double e_;
		double sece_;

		friend class TransverseMercator;
	};

	class TransverseMercator {
	public:
		TransverseMercator(Ellipsoid ellipsoid=Ellipsoid(), double lon0=0,
			double lat0 = 0, double k0 = 1, double fe = 500000, double fn = 0);
		TransverseMercator(const TransverseMercator& proj);
		TransverseMercator& operator=(const TransverseMercator& proj);

		void GeodeticToProjection(double lat, double lon, double& east, double& north);
		void ProjectionToGeodetic(double east, double north, double& lat, double& lon);

	protected:
		Ellipsoid ellipsoid_;
		double Latitude_Of_Origin;
		double Central_Meridian;
		double Scale_Factor;
		double False_Easting;
		double False_Northing;
	};

}


#endif
