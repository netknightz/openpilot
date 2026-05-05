from cereal import log
from openpilot.system.hardware.base import HardwareBase

NetworkType = log.DeviceState.NetworkType


class Pc(HardwareBase):
  def get_device_type(self):
    return "tizi"

  def get_network_type(self):
    return NetworkType.wifi
