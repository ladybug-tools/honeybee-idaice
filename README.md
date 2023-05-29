# Honeybee -> IDA-ICE

A honeybee extension to convert HBJSON files to [IDA ICE](https://www.equa.se/en/ida-ice) `idm` files.

![image](https://github.com/ladybug-tools/honeybee-idaice/assets/2915573/12840ca3-a94f-4c68-9251-e8eb393bc048)

Two comments on how the model should be prepared.

1. IDA ICE expects the model to be exported at the finish line of the one - and not the center of the wall.
1. IDA ICE intersects the faces automatically. Do not intersect the faces in HBJSON files.

The exporter only exports the geometry. None of the energy or Radiance properties are exported.
