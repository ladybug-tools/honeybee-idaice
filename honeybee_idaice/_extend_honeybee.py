from .writer import model_to_idm
from honeybee.model import Model

Model.to_idm = model_to_idm
