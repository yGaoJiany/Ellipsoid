
#define DLL_ELLIPSOID_EXPORTS

#include <Eigen/Eigen>
#include "Ellipsoid.h"


namespace coti {
	Ellipsoid::Ellipsoid(double a /*= 6378137.0*/, double invf /*= 298.257223563*/) {
		a_ = a;
		invf_ = invf;
		f_ = 1.0 / invf_;
		b_ = a_ * (1 - f_);
		e_ = sqrt(a_ * a_ - b_ * b_) / a_;
		sece_ = sqrt(a_ * a_ - b_ * b_) / b_;
	}

	Ellipsoid::Ellipsoid(const Ellipsoid& ellipsoid) {
		this->a_ = ellipsoid.a_;
		this->invf_ = ellipsoid.invf_;
		this->f_ = ellipsoid.f_;
		this->b_ = ellipsoid.b_;
		this->e_ = ellipsoid.e_;
		this->sece_ = ellipsoid.sece_;
	}

	Ellipsoid& Ellipsoid::operator=(const Ellipsoid& ellipsoid) {
		if (this == &ellipsoid) {
			return *this;
		}
		this->a_ = ellipsoid.a_;
		this->invf_ = ellipsoid.invf_;
		this->f_ = ellipsoid.f_;
		this->b_ = ellipsoid.b_;
		this->e_ = ellipsoid.e_;
		this->sece_ = ellipsoid.sece_;

		return *this;
	}

	double radian2degree(double rad) {
		return 180 / acos(-1) * rad;
	}

	double degree2radian(double deg) {
		return acos(-1) / 180 * deg;
	}

	double Ellipsoid::Z_XbarPsi(double xBar, double psi) {
		double tmp;
		tmp = sqrt(b_ * b_ / (1.0 - e_ * e_ * cos(psi) * cos(psi)));
		tmp *= (a_ * a_ / (b_ * b_) - 1.0) * sin(psi);
		return(xBar * a_ * a_ / (b_ * b_) * tan(psi) - tmp);
	}

	bool Ellipsoid::GeocentricToGeodetic(
		double X, double Y, double Z,
		double& lat, double& lon, double& hei) {

		double  delta = 0.00001;
		double xBar, psi, re, z1, dz;

		lon = atan2(Y, X);

		xBar = sqrt(X * X + Y * Y);

		psi = atan2(Z, xBar);

		int num = 0;

		do {
			z1 = Z_XbarPsi(xBar, psi);
			dz = Z_XbarPsi(xBar, psi + delta) - z1;

			psi += (Z - z1) / dz * delta;
			num++;
			if (num > 1000)
				break;
		} while (fabs(Z - z1) > fabs(dz));


		lat = atan(tan(psi) * a_ * a_ / (b_ * b_));

		re = sqrt(b_ * b_ / (1.0 - e_ * e_ * cos(psi) * cos(psi)));

		hei = sqrt((Z - re * sin(psi)) * (Z - re * sin(psi))
			+ (xBar - re * cos(psi)) * (xBar - re * cos(psi)));

		lat = radian2degree(lat);
		lon = radian2degree(lon);

		if ((X * X + Y * Y + Z * Z) < re * re)
			hei *= -1.0;
		if (num > 1000)
			return false;
		else
			return true;
	}

	void Ellipsoid::GeocentricToGeodeticOzone(
		double X, double Y, double Z,
		double& lat, double& lon, double& hei) {
		int is_north = 1;
		double phi = 0;
		double p = sqrt(X * X + Y * Y);
		if (Z != 0) {
			is_north = (Z > 0) ? 1 : -1;
			Z *= is_north;
			double c2 = e_ * e_ * a_ * a_;
			double N = (a_ * p - c2) / (2 * b_ * Z);
			double S = (a_ * p + c2) / (2 * b_ * Z);

			double V = (4 * N * S + 1) / 3;
			double W = S * S - N * N;
			double D = sqrt(W * W + V * V * V);
			double I = pow(D + W, 1.0 / 3.0) - pow(D - W, 1.0 / 3.0);
			double J = sqrt(4 * N * N + 2 * I);
			double K = sqrt(I * I + 1);
			double T = 2 * N + J;
			double G = T * T - 4 * (I - K);
			double U = (sqrt(G) + T) / 2.0;
			phi = atan2(2 * a_ * U, b_ * (U * U - 1));
		}

		Z *= is_north;

		lat = phi * is_north;
		lon = atan2(Y, X);
		hei = p / cos(lat) - a_ * a_ / sqrt(a_ * a_ * cos(lat) * cos(lat) + b_ * b_ * sin(lat) * sin(lat));

		lat = radian2degree(lat);
		lon = radian2degree(lon);
	}

	void Ellipsoid::GeodeticToGeocentric(
		double lat, double lon, double hei,
		double& X, double& Y, double& Z) {

		lat = degree2radian(lat);
		lon = degree2radian(lon);

		double re, psi, xBar;

		psi = atan2(b_ * b_ * tan(lat), a_ * a_);

		re = sqrt(b_ * b_ / (1.0 - e_ * e_ * cos(psi) * cos(psi)));

		xBar = re * cos(psi) + hei * cos(lat);

		X = xBar * cos(lon);
		Y = xBar * sin(lon);
		Z = re * sin(psi) + hei * sin(lat);
	}

	void Ellipsoid::ECEFToENU(double X, double Y, double Z, double lat0, double lon0, double h0, double& x, double& y, double& z) {
		double xo, yo, zo;

		GeodeticToGeocentric(lat0, lon0, h0, xo, yo, zo);

		Eigen::Vector3d vect(X - xo, Y - yo, Z - zo);

		lat0 = degree2radian(lat0);
		lon0 = degree2radian(lon0);

		double sinlat = sin(lat0);
		double coslat = cos(lat0);
		double sinlon = sin(lon0);
		double coslon = cos(lon0);

		Eigen::Matrix<double, 3, 3> Rotation;
		Rotation(0, 0) = -sinlon; Rotation(0, 1) = coslon; Rotation(0, 2) = 0;
		Rotation(1, 0) = -sinlat * coslon; Rotation(1, 1) = -sinlat * sinlon; Rotation(1, 2) = coslat;
		Rotation(2, 0) = coslat * coslon; Rotation(2, 1) = coslat * sinlon; Rotation(2, 2) = sinlat;

		Eigen::Vector3d xyz = Rotation * vect;

		x = xyz[0];
		y = xyz[1];
		z = xyz[2];
	}

	void Ellipsoid::ENUToECEF(double x, double y, double z, double lat0, double lon0, double h0, double& X, double& Y, double& Z) {
		Eigen::Vector3d vect(x, y, z);

		double xo, yo, zo;
		GeodeticToGeocentric(lat0, lon0, h0, xo, yo, zo);

		lat0 = degree2radian(lat0);
		lon0 = degree2radian(lon0);
		double sinlat = sin(lat0);
		double coslat = cos(lat0);
		double sinlon = sin(lon0);
		double coslon = cos(lon0);

		Eigen::Matrix<double, 3, 3> Rotation;
		Rotation(0, 0) = -sinlon; Rotation(0, 1) = -sinlat * coslon;  Rotation(0, 2) = coslat * coslon;
		Rotation(1, 0) = coslon; Rotation(1, 1) = -sinlat * sinlon; Rotation(1, 2) = coslat * sinlon;
		Rotation(2, 0) = 0;  Rotation(2, 1) = coslat; Rotation(2, 2) = sinlat;

		Eigen::Vector3d XYZ = Rotation * vect;

		X = XYZ[0] + xo;
		Y = XYZ[1] + yo;
		Z = XYZ[2] + zo;
	}

	TransverseMercator::TransverseMercator(Ellipsoid ellipsoid, double lon0, double lat0 /*= 0*/, double k0 /*= 1*/, double fe /*= 500000*/, double fn /*= 0*/)
	{
		ellipsoid_ = ellipsoid;
		Central_Meridian = lon0;
		Latitude_Of_Origin = lat0;
		Scale_Factor = k0;
		False_Easting = fe;
		False_Northing = fn;
	}

	TransverseMercator::TransverseMercator(const TransverseMercator& proj) {
		this->ellipsoid_ = proj.ellipsoid_;
		this->Latitude_Of_Origin = proj.Latitude_Of_Origin;
		this->Central_Meridian = proj.Central_Meridian;
		this->Scale_Factor = proj.Scale_Factor;
		this->False_Easting = proj.False_Easting;
		this->False_Northing = proj.False_Northing;
	}

	TransverseMercator& TransverseMercator::operator=(const TransverseMercator& proj) {
		if (this == &proj) {
			return *this;
		}
		this->ellipsoid_ = proj.ellipsoid_;
		this->Latitude_Of_Origin = proj.Latitude_Of_Origin;
		this->Central_Meridian = proj.Central_Meridian;
		this->Scale_Factor = proj.Scale_Factor;
		this->False_Easting = proj.False_Easting;
		this->False_Northing = proj.False_Northing;

		return *this;
	}

	void TransverseMercator::GeodeticToProjection(double lat, double lon, double& east, double& north) {
		lat = degree2radian(lat);
		lon = degree2radian(lon);

		double lat0 = degree2radian(Latitude_Of_Origin);
		double lon0 = degree2radian(Central_Meridian);

		//Calculate Then the meridional arc distance from equator to the projection origin (M0)
		double e2 = ellipsoid_.e_ * ellipsoid_.e_;
		double e4 = e2 * e2;
		double e6 = e4 * e2;

		double M0 = ellipsoid_.a_ * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lat0 -
			(3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * sin(2 * lat0) +
			(15 * e4 / 256 + 45 * e6 / 1024) * sin(4 * lat0) - (35 * e6 / 3072) * sin(6 * lat0));

		//calculate T C A v M
		double cos_lat = cos(lat);
		double sin_lat = sin(lat);
		double tan_lat = tan(lat);

		double T = tan_lat * tan_lat;
		double C = e2 * cos_lat * cos_lat / (1 - e2);
		double A = (lon - lon0) * cos_lat;
		double v = ellipsoid_.a_ / sqrt(1 - e2 * sin_lat * sin_lat);
		double M = ellipsoid_.a_ * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lat -
			(3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * sin(2 * lat) +
			(15 * e4 / 256 + 45 * e6 / 1024) * sin(4 * lat) - (35 * e6 / 3072) * sin(6 * lat));

		double A2 = A * A;
		double A3 = A * A * A;

		east = False_Easting + Scale_Factor * v * (A + (1 - T + C) * A3 / 6 +
			(5 - 18 * T + T * T + 72 * C - 58 * e2) * A2 * A3 / 120);
		north = False_Northing + Scale_Factor * (M - M0 + v * tan_lat * (A2 / 2 + (5 - T + 9 * C + 4 * C * C) * A2 * A2 / 24 +
			(61 - 58 * T + T * T + 600 * C - 330 * e2) * A3 * A3 / 720));
	}

	void TransverseMercator::ProjectionToGeodetic(double east, double north, double& lat, double& lon) {
		double lat0 = degree2radian(Latitude_Of_Origin);
		double lon0 = degree2radian(Central_Meridian);

		double e2 = ellipsoid_.e_ * ellipsoid_.e_;
		double e4 = e2 * e2;
		double e6 = e4 * e2;

		double	M0 = ellipsoid_.a_ * ((1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256) * lat0 -
			(3 * e2 / 8 + 3 * e4 / 32 + 45 * e6 / 1024) * sin(2 * lat0) +
			(15 * e4 / 256 + 45 * e6 / 1024) * sin(4 * lat0) - (35 * e6 / 3072) * sin(6 * lat0));

		// calculate e1 u1 M1
		double temp_e = sqrt(1 - ellipsoid_.e_ * ellipsoid_.e_);
		double e1 = (1 - temp_e) / (1 + temp_e);
		double M1 = M0 + (north - False_Northing) / Scale_Factor;
		double u1 = M1 / (ellipsoid_.a_ * (1 - e2 / 4 - 3 * e4 / 64 - 5 * e6 / 256));

		//calculate lat1
		double	e1_2 = e1 * e1;
		double	lat1 = u1 + (3 * e1 / 2 - 27 * e1_2 * e1 / 32) * sin(2 * u1) +
			(21 * e1_2 / 16 - 55 * e1_2 * e1_2 / 32) * sin(4 * u1) +
			(151 * e1_2 * e1 / 96) * sin(6 * u1) + (1097 * e1_2 * e1_2 / 512) * sin(8 * u1);

		double	temp = sqrt(1 - e2 * sin(lat1) * sin(lat1));
		double	v1 = ellipsoid_.a_ / temp;
		double	p1 = ellipsoid_.a_ * (1 - e2) / (temp * temp * temp);
		double	T1 = tan(lat1) * tan(lat1);

		double	C1 = ellipsoid_.sece_ * cos(lat1);
		C1 = C1 * C1;

		double D = (east - False_Easting) / (v1 * Scale_Factor);

		//calculate lat, lon
		double	D2 = D * D;
		double	D3 = D2 * D;

		double sece2 = ellipsoid_.sece_ * ellipsoid_.sece_;

		lat = lat1 - (v1 * tan(lat1) / p1) * (D2 / 2 - (5 + 3 * T1 + 10 * C1 - 4 * C1 * C1 - 9 * sece2) * D2 * D2 / 24 + (
			61 + 90 * T1 + 298 * C1 + 45 * T1 * T1 - 252 * sece2 - 3 * C1 * C1) * D3 * D3 / 720);
		lon = lon0 + (D - (1 + 2 * T1 + C1) * D3 / 6 + (5 - 2 * C1 + 28 * T1 - 3 * C1 * C1 + 8 * sece2 + 24 * T1 * T1) *
			D2 * D3 / 120) / cos(lat1);

		Eigen::Vector2d BL;
		lat = radian2degree(lat);
		lon = radian2degree(lon);
	}
}
