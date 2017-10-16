""" test gb_box. """

import box

INFO = {'grid_file_name': '/Users/CnWang/Documents/gb_roms/grd/GlacierBay_lr_grd.nc',
        'river_raw_file_name': '/Users/CnWang/Documents/gb_roms/gb_discharge.nc',
        'river_file_name': '/Users/CnWang/git/gb_box/data/gb_box_rivers.nc',
        'wind_file_name': '/Users/CnWang/git/gb_box/data/juneau.csv',
        'sp_raw_file_name': '/Users/CnWang/Documents/gb_roms/CTDS7513',
        'sp_file_name': '/Users/CnWang/git/gb_box/data/gb_box_sp.nc',
        'sl_rivers': 'l',
        'sl_sp': 'l'}

box_info = box.Box(INFO)
box_info()
rivers = box.BoxRivers(box_info)
rivers()
cdo = box.BoxCDO(box_info)
cdo()
spraw = box.BoxSp(box_info)
spraw()
