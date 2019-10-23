"""
Imports all the functions in the various files.
Doing this explicitly for clarity
"""

from .array_help import (
    find_nearest,
    find_nearest_index,
    copy_range,
    copy_range_array,
    copy_range_yarray,
    sort_by_list,
    find_nearest_tolerance,
    sort_array_column,
    split_columns,
)

from .background import subtract_background, find_background_lp_filter

from .peaks import (
    find_peaks,
    find_peaks_pu,
    lorentzian,
    gaussian,
    voigt,
    lorentzian_d,
    lorentzian_dd,
    gaussian_d,
    gaussian_dd,
)

from .fitting import (
    make_region_array,
    make_regions,
    fit_data,
    fit_data_bg,
    output_results,
    batch_fit_single_peak,
    line_fit,
    exponential_fit,
    exponential_fit_offset,
    poly_fit,
    fit_peaks,
    split_and_fit,
    fit_peak_table,
    build_model,
    set_parameters,
    build_model_d,
    build_model_dd,
    peak_table,
)

from .convert import (
    wl2wn,
    wl2rwn,
    wn2wl,
    rwn2wl,
    rwn2wn,
    absorption,
    nm2ev,
    nm2ev_xy,
    nm2ev_xyz,
)
from .file_io import (
    assure_path_exists,
    write2col,
    getxy,
    clean_file,
    path,
    load_folder,
    quick_load_xy,
    make_fname,
    rpath,
    read_json,
    write_json,
)
from .normalize import (
    normalize,
    normalize_msc,
    normalize_pq,
    normalize_2pt,
    normalize_fs,
)
from .filters import smooth_data, butter_lp_filter
from .calibrate import (
    calibrate,
    find_laser_wavelength,
    find_best_offset,
    find_best_offset2,
    calibrate_x_data2,
    neon_peaks,
    calibrate_neon_wavenumber,
)
from .read_files import read_cary, read_craic, read_nicolet, read_horiba, read_renishaw
from .misc import remove_absorption_jumps, generate_spectrum, activity_to_intensity

from .plot import plot_peak_fit, plot_components, plot_background, plot_fit

from .filters import resample

# Need to actually make some of these work
# from filters import smooth_data, butter_lp_filter, butter_lowpass_filter, \
#                     butter_lowpass, wicker, savgol, butterworth_bandpass
