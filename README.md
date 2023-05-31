# Honeybee -> IDA-ICE

A honeybee extension to convert HBJSON files to [IDA ICE](https://www.equa.se/en/ida-ice) `idm` files.

![Revit Sample Model](https://github.com/ladybug-tools/honeybee-idaice/assets/2915573/97ce39b6-8f45-4dfc-b2f6-152d457a9c82)

Two comments on how the model should be prepared.

1. IDA ICE expects the model to be exported at the finish line of the one instead of the center of the wall.
1. IDA ICE intersects the faces automatically. Do not intersect the faces in HBJSON files.
If the model is already intersected, use Pollination Rhino's `PO_RebuildRooms` command to merge the faces together.

The exporter only exports the geometry. None of the energy or Radiance properties are exported.
