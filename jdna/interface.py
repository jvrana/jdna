from abc import ABC
from abc import abstractmethod


class ClassInterface(ABC):
    def __init__(self, cls):
        self._cls = cls

    @abstractmethod
    def instance(self, instance):
        pass


class Interface:
    def __init__(self, instance, class_interface):
        self._inst = instance
        self._class_interface = class_interface
