# coding=utf-8
# import all of the modules for writing geometry to IDAICE
from honeybee.properties import ModelProperties

from .properties.model import ModelIDAICEProperties

# set a hidden ies attribute on each core geometry Property class to None
# define methods to produce ies property instances on each Property instance
ModelProperties._idaice = None


def model_idaice_properties(self):
    if self._idaice is None:
        self._idaice = ModelIDAICEProperties(self.host)
    return self._idaice


# add IDAICE property methods to the Properties classes
ModelProperties.idaice = property(model_idaice_properties)
