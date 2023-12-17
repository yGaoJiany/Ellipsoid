
# Coordinate Conversions and Transformation
There are many coordinate systems in satellite photogrammetry, such as Geodetic Coordinate System, Projected Coordinate System, Scene East-North-Up (Scene ENU) Coordinate System and so on. We only give a brief introduction to some widely used coordinate systems here.

This repo implements the following coordinate conversions:
- Conversion between **Geographic** and **Geocentric**
- Conversion between **Geodetic** and **Projection** (Transverse Mercator Projection)
- Conversion between **ECEF** and **ENU**

The test data are generated using the [geodetic-calculator](https://www.midpointgeo.net/geodetic-calculator/).

[Eigen](https://eigen.tuxfamily.org/index.php?title=Main_Page) is required to compile the code.

## More references
[1] Coordinate Conversions and Transformations including Formulas. https://www.iogp.org/wp-content/uploads/2019/09/373-07-02.pdf

[2] For more Projected Coordinate Systems, please refer to https://developers.arcgis.com/javascript/3/jshelp/pcs.htm

[3] For more Geographic Coordinate Systems, please refer to https://developers.arcgis.com/javascript/3/jshelp/gcs.htm