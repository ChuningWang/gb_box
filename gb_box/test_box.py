""" test gb_box. """

import box

INFO = {'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
        'river_raw_file_name': '/Users/CnWang/Documents/gb_roms/gb_discharge.nc',
        'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers.nc',
        'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
        'sp_raw_file_name': '/Users/CnWang/Documents/gb_roms/CTDS7513',
        'sp_soda_dir_name': '/Users/CnWang/Documents/gb_roms/CTDS7513',
        'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc',
        'sl_rivers': 'l',
        'sl_sp': 'l'}

INFO = {'grid_file_name': '/Users/chuning/GlacierBay_lr_grd.nc',
        'river_raw_file_name': '/glade/p/work/chuning/data/gb_discharge.nc',
        'river_file_name': '/glade/u/home/chuning/git/gb_box/data/gb_box_rivers.nc',
        'wind_file_name': '/glade/u/home/chuning/git/gb_box/data/juneau.csv',
        'sp_soda_dir_name': '/Volumes/P1/Data/SODA/SODA_3.3.1/',
        'sp_file_name': '/Users/chuning/git/gb_box/data/gb_box_sp.nc',
        'sl_rivers': 's',
        'sl_sp': 's'}

box_info = box.Box(INFO)
box_info()
# rivers = box.BoxRivers(box_info)
# rivers()
# cdo = box.BoxCDO(box_info)
# cdo()
spraw = box.BoxSp(box_info)
spraw()
