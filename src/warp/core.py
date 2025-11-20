from . import cascades
from . import dace


class Star:
    def __init__(self, name, instrument=None, load_tess=False, max_erv=30, adjust_means=True, do_secular_corr=True):
        self.name = name
        self.rv_data = dace.download_points(name)
