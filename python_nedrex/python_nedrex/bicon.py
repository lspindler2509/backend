from pathlib import Path as _Path
from typing import IO as _IO
from typing import Any as _Any
from typing import Dict as _Dict
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import check_response as _check_response
from python_nedrex.common import download_file as _download_file
from python_nedrex.common import http as _http


def bicon_request(
    expression_file: _IO[str],
    lg_min: int = 10,
    lg_max: int = 15,
    network: str = "DEFAULT",
) -> str:
    files = {"expression_file": expression_file}
    data = {"lg_min": lg_min, "lg_max": lg_max, "network": network}

    url = f"{_config.url_base}/bicon/submit"
    resp = _http.post(url, data=data, files=files, headers={"x-api-key": _config.api_key})
    result: str = _check_response(resp)
    return result


def check_bicon_status(uid: str) -> _Dict[str, _Any]:
    url = f"{_config.url_base}/bicon/status"
    resp = _http.get(url, params={"uid": uid}, headers={"x-api-key": _config.api_key})
    result: _Dict[str, _Any] = _check_response(resp)
    return result


def download_bicon_data(uid: str, target: _Optional[str] = None) -> None:
    if target is None:
        target = str(_Path(f"{uid}.zip").resolve())

    url = f"{_config.url_base}/bicon/download?uid={uid}"

    _download_file(url, target)
