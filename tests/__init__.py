import os
import sys

from pathlib import Path

sys.path.append(os.path.abspath(
    os.path.join(os.path.realpath(__file__), '../../')
))

TEST_DATA_DIR = f"{Path(__file__).parent.resolve()}/test_data"