"""Top-level package for python-nedrex."""

__author__ = """David James Skelton"""
__email__ = "james.skelton@newcastle.ac.uk"
__version__ = "0.1.1"

from typing import Optional

from attrs import define


@define
class _Config:
    _url_base: Optional[str] = None
    _url_vpd: Optional[str] = None
    _api_key: Optional[str] = None

    @property
    def url_base(self) -> Optional[str]:
        """Returns the API URL base stored on the _Config instance"""
        return self._url_base

    @property
    def url_vpd(self) -> Optional[str]:
        """Returns the VPD URL base stored on the _Config instance"""
        return self._url_vpd

    @property
    def api_key(self) -> Optional[str]:
        """Returns the API key stored on the _Config instance"""
        return self._api_key

    def set_url_base(self, url_base: str) -> None:
        """Sets the URL base for the API in the configuration"""
        self._url_base = url_base.rstrip("/")

    def set_url_vpd(self, url_vpd: str) -> None:
        """Sets the URL base for the VPD in the configuration"""
        self._url_vpd = url_vpd.rstrip("/")

    def set_api_key(self, key: str) -> None:
        """Sets the API key in the configuration"""
        self._api_key = key


config: _Config = _Config()
