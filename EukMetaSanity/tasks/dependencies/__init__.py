"""
Populate dependencies for easy loading
"""
from EukMetaSanity.tasks.utils.load_modules import get_modules

# # TODO: Add `expects` to docstrings for all dependencies
dependencies = get_modules(__file__)
