
#include<iostream>
#include"Ellipsoid.h"

// The test data is calculated from https://www.midpointgeo.net/ Geodetic Calculator
void EXPECT_LT(double x, double y){
	if(x > y)
		std::cout << "Failed! " << x << " is lager than "<< y <<std::endl;
	else
		std::cout << "Pass the test." <<std::endl;
}

void test_ELLIPSOID_GEOCENTRIC2GEODETIC() {
	std::cout << "test_ELLIPSOID_GEOCENTRIC2GEODETIC" << std::endl;

	double X = -3010904.36;
	double Y = 4636386.14;
	double Z = 3170373.74;

	coti::Ellipsoid wgs84;
	double lat, lon, hei;
	wgs84.GeocentricToGeodetic(X, Y, Z, lat, lon, hei);
	
	EXPECT_LT(fabs(lat - 30), 1e-7);
	EXPECT_LT(fabs(lat - 30), 1e-7);
	EXPECT_LT(fabs(hei), 1e-2);
}

void test_ELLIPSOID_GEOCENTRIC2GEODETICOZONE() {
	std::cout << "test_ELLIPSOID_GEOCENTRIC2GEODETICOZONE" << std::endl;

	double X = -3010904.36;
	double Y = 4636386.14;
	double Z = 3170373.74;

	coti::Ellipsoid wgs84;
	double lat, lon, hei;
	wgs84.GeocentricToGeodeticOzone(X, Y, Z, lat, lon, hei);

	EXPECT_LT(fabs(lat - 30), 1e-7);
	EXPECT_LT(fabs(lon - 123), 1e-7);
	EXPECT_LT(fabs(hei), 1e-2);
}

void test_ELLIPSOI_GEODETIC2GEOCENTRIC() {
	std::cout << "test_ELLIPSOI_GEODETIC2GEOCENTRIC" <<std::endl;
	double lat = 30;
	double lon = 123;
	double hei = 0;

	coti::Ellipsoid wgs84;
	double X, Y, Z;
	wgs84.GeodeticToGeocentric(lat, lon, hei, X, Y, Z);

	EXPECT_LT(fabs(X - (-3010904.36)), 1e-2);
	EXPECT_LT(fabs(Y - 4636386.14), 1e-2);
	EXPECT_LT(fabs(Z - 3170373.74), 1e-2);
}

void test_ELLIPSOID_ECEF2ENU() {
	std::cout <<"test_ELLIPSOID_ECEF2ENU"<<std::endl;

	double X = 660930;
	double Y = -4701424;
	double Z = 4246579;

	coti::Ellipsoid wgs84;
	double E, N, U;
	double lat0, lon0, hei0;
	lat0 = 42; lon0 = -82; hei0 = 200;
	wgs84.ECEFToENU(X, Y, Z, lat0, lon0, hei0, E, N, U);

	EXPECT_LT(fabs(E - 186.12), 1e-2);
	EXPECT_LT(fabs(N - 286.56), 1e-2);
	EXPECT_LT(fabs(U - 939.10), 1e-2);
}

void test_ELLIPSOIDTEST_ENU2ECEF() {
	std::cout <<"test_ELLIPSOIDTEST_ENU2ECEF"<<std::endl;

	double E = 186.12;
	double N = 286.56;
	double U = 939.10;

	coti::Ellipsoid wgs84;
	double X, Y, Z;
	wgs84.ENUToECEF(E, N, U, 42, -82, 200, X, Y, Z);

	EXPECT_LT(fabs(X - 660930), 1e-2);
	EXPECT_LT(fabs(Y + 4701424), 1e-2);
	EXPECT_LT(fabs(Z - 4246579), 1e-2);
}

void test_TRANSVERSEMERCATOR_GEODETIC2PROJECTION() {
	std::cout <<"test_TRANSVERSEMERCATOR_GEODETIC2PROJECTION"<<std::endl;

	double lat = 30;
	double lon = 120;

	coti::Ellipsoid wgs84;
	coti::TransverseMercator proj(wgs84, 123, 0, 0.9996);
	double E, N;
	proj.GeodeticToProjection(lat, lon, E, N);

	EXPECT_LT(fabs(E - 210590.35), 1e-2);
	EXPECT_LT(fabs(N - 3322575.90), 1e-2);
}

void test_TRANSVERSEMERCATOR_PROJECTION2GEODETIC() {
	std::cout <<"test_TRANSVERSEMERCATOR_PROJECTION2GEODETIC"<<std::endl;
	
	double X = 210590.35;
	double Y = 3322575.90;

	coti::Ellipsoid wgs84;
	coti::TransverseMercator proj(wgs84, 123, 0, 0.9996);
	double lat, lon;
	proj.ProjectionToGeodetic(X, Y, lat, lon);

	EXPECT_LT(fabs(lat - 30), 1e-7);
	EXPECT_LT(fabs(lon - 120), 1e-7);
}

int main(){
	test_ELLIPSOID_GEOCENTRIC2GEODETIC();
	test_ELLIPSOID_GEOCENTRIC2GEODETICOZONE();
	test_ELLIPSOI_GEODETIC2GEOCENTRIC();
	test_ELLIPSOID_ECEF2ENU();
	test_ELLIPSOIDTEST_ENU2ECEF();
	test_TRANSVERSEMERCATOR_GEODETIC2PROJECTION();
	test_TRANSVERSEMERCATOR_PROJECTION2GEODETIC();
}