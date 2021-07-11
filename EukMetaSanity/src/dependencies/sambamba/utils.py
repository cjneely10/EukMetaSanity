# Detect binary file format for SAM/BAM
# Seen on https://stackoverflow.com/questions/898669/how-can-i-detect-if-a-file-is-binary-non-text-in-python
from pathlib import Path
from typing import Union

textchars = bytearray({7, 8, 9, 10, 12, 13, 27} | set(range(0x20, 0x100)) - {0x7f})


def is_binary_string(_bytes: Union[str, bytes]) -> bool:
    return bool(_bytes.translate(None, textchars))


def is_sam(file: Path) -> bool:
    return not is_binary_string(open(file, "rb").read(1024))
