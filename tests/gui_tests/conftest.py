from pathlib import Path
import pytest

@pytest.fixture
def ROOT_DIR() -> Path:
    return Path(__file__).parent.parent.parent 

@pytest.fixture
def TEST_DATA_DIR() -> Path:
    return Path(__file__).parent / "test_data"
