# -*- coding: utf-8 -*-
from cellrank.utils.models import Model
from typing import Union, Mapping


_model_type = Union[Model, Mapping[str, Mapping[str, Model]]]
