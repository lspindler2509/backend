import os as _os
from typing import Optional as _Optional

from python_nedrex import config as _config
from python_nedrex.common import http as _http
from python_nedrex.decorators import check_url_vpd as _check_url_vpd


@_check_url_vpd
def get_vpd(disorder: str, number_of_patients: int, out_dir: str) -> _Optional[str]:
    """
    Downloads a .zip archive with the requested virtual patient data to the given directory.

        Parameters:
            disorder (str): The disorder mondoID (e.g mondo.0000090) for which the virtual patient should be retrieved.
            number_of_patients (int): The number of simulated patients in the dataset. Can be 1, 10 or 100.
            out_dir (str): The absolute path of a directory where the virtual patient data should be stored.

        Returns:
            archive (str): Absolute path of the downloaded zip archive or None if the requested resource does not exist.
    """
    archive_name: str = f"{disorder}_1000GP_{number_of_patients}VPSim.zip"
    url: str = f"{_config.url_vpd}/vpd/{disorder}/{archive_name}"
    archive: str = _os.path.join(out_dir, archive_name)

    data = _http.get(url)
    if data.status_code != 200:
        return None

    with open(archive, "wb") as arch:
        arch.write(data.content)
    return archive
