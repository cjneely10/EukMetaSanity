#!/usr/bin/env python3
import os
from EukMetaSanity.utils.arg_parse import ArgParse


if __name__ == "__main__":
    _ap = ArgParse(
        (
            ((), {}),
            ((), {}),
            ((), {}),
            ((), {}),
        ),
        description="Convert BED to GFF3 format"
    )
