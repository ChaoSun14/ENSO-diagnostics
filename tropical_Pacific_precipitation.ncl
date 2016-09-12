load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/shea_util.ncl"
begin
    
    OUTPUT_TYPE = "png"
    OUTDIR      = "/home/reproduce/sunc/enso/plot/"
    CESMDIR     = "/home/reproduce/sunc/enso/CESM1-FASTCHEM_picontrol_atm/"
    FGOALSDIR   = "/home/reproduce/sunc/enso/FGOALS-g2_picontrol_atm/"
    GPCP        = "/home/reproduce/sunc/enso/obs_data/"

    wks_type = OUTPUT_TYPE

    wks_pre = gsn_open_wks(wks_type,OUTDIR+"tropical_pre")
    gsn_define_colormap(wks_pre,"amwg")


    data = addfile(GPCP+"precip.mon.mean.nc","r") ;197901-201503
    precip_obs = data->precip(0:36*12-1,:,:)    ;197901-201412 mm/day
    lat_obs = data->lat
    lon_obs = data->lon
    time_obs = data->time(0:36*12-1)
    ;printVarSummary(precip_obs)

    pi=4.*atan(1.0)
    rad=(pi/180.)
    coswgt_obs=cos(rad*lat_obs)
    coswgt_obs!0 = "lat" 
    coswgt_obs&lat= lat_obs 

    ;d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    ;basemap = d->LSMASK
    ;lsm = landsea_mask(basemap,lat_obs,lon_obs)
    ;precip_obs = mask(precip_obs,conform(precip_obs,lsm,(/1,2/)).ge.1,False)
    ;delete(lsm)

    precip_obs_mon = dim_avg_n_Wrap(precip_obs, 0)
    precip_obs_mon_ave = wgt_areaave_Wrap(precip_obs_mon({-5:5},{120:290}),coswgt_obs({-5:5}),1.0,1)
    ;print(precip_obs_mon_ave)


    data = addfile(CESMDIR+"pr_Amon_CESM1-FASTCHEM_piControl_r1i1p1_007001-029112.nc","r") 
    precip_cesm = data->pr(0:200*12-1,:,:)    
    lat_cesm = data->lat
    lon_cesm = data->lon
    time_cesm = data->time(0:200*12-1)
    printVarSummary(precip_cesm)

    coswgt_cesm=cos(rad*lat_cesm)
    coswgt_cesm!0 = "lat" 
    coswgt_cesm&lat= lat_cesm 
    precip_cesm=(/precip_cesm*24*60*60/)   ;kg m-2 s-1 to mm/day
    precip_cesm@units = "mm/day"

    ;d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    ;basemap = d->LSMASK
    ;lsm = landsea_mask(basemap,lat_cesm,lon_cesm)
    ;precip_cesm = mask(precip_cesm,conform(precip_cesm,lsm,(/1,2/)).ge.1,False)
    ;delete(lsm)

    precip_cesm_mon = dim_avg_n_Wrap(precip_cesm, 0)
    precip_cesm_mon_ave = wgt_areaave_Wrap(precip_cesm_mon({-5:5},{120:290}),coswgt_cesm({-5:5}),1.0,1)
    ;print(precip_cesm_mon_ave)
    precip_cesm_mon_obsgrid = linint2_Wrap(lon_cesm, lat_cesm, precip_cesm_mon, True, lon_obs, lat_obs, 0)
    ;printVarSummary(precip_cesm_mon_obsgrid)
    precip_cesm_mon_rmse = dim_rmsd_n_Wrap(precip_cesm_mon_obsgrid({-5:5},{120:290}),precip_obs_mon({-5:5},{120:290}),(/0,1/))
    ;print(precip_cesm_mon_rmse)

    data = addfile(FGOALSDIR+"t0225_ocn_S.gamil.h0.PRECT.060001-079912.nc","r") 
    precip_fgoals = data->PRECT(0:200*12-1,:,:)    
    lat_fgoals = data->lat
    lon_fgoals = data->lon
    time_fgoals = data->time(0:200*12-1)
    ;printVarSummary(precip_fgoals)

    coswgt_fgoals=cos(rad*lat_fgoals)
    coswgt_fgoals!0 = "lat" 
    coswgt_fgoals&lat= lat_fgoals 
    precip_fgoals=(/precip_fgoals*24*60*60*1000/)   ;m/s to mm/day
    precip_fgoals@units = "mm/day"

    ;d = addfile("$NCARG_ROOT/lib/ncarg/data/cdf/landsea.nc","r")   
    ;basemap = d->LSMASK
    ;lsm = landsea_mask(basemap,lat_fgoals,lon_fgoals)
    ;precip_fgoals = mask(precip_fgoals,conform(precip_fgoals,lsm,(/1,2/)).ge.1,False)
    ;delete(lsm)

    precip_fgoals_mon = dim_avg_n_Wrap(precip_fgoals, 0)
    precip_fgoals_mon_ave = wgt_areaave_Wrap(precip_fgoals_mon({-5:5},{120:290}),coswgt_fgoals({-5:5}),1.0,1)
    ;print(precip_fgoals_mon_ave)
    precip_fgoals_mon_obsgrid = linint2_Wrap(lon_fgoals, lat_fgoals, precip_fgoals_mon, True, lon_obs, lat_obs, 0)
    ;printVarSummary(precip_fgoals_mon_obsgrid)
    precip_fgoals_mon_rmse = dim_rmsd_n_Wrap(precip_fgoals_mon_obsgrid({-5:5},{120:290}),precip_obs_mon({-5:5},{120:290}),(/0,1/))
    ;print(precip_fgoals_mon_rmse)


    plot1 = new(3,graphic)
    res1=True
    res1@gsnDraw        = False            ; don't draw yet
    res1@gsnFrame       = False
    ;res1@mpCenterLonF   = 180.
    res1@vpHeightF      = 0.3
    res1@vpWidthF       = 0.7
    res1@gsnStringFontHeightF      = 0.022
    res1@gsnLeftStringFontHeightF  = 0.022
    res1@gsnRightStringFontHeightF = 0.022
    res1@tmXBLabelFontHeightF  = 0.022
    res1@tmLabelAutoStride     = True
    res1@tmYLLabelFontHeightF  = 0.022
    res1@cnFillOn  = True
    res1@cnLinesOn = False
    res1@lbLabelBarOn = False
    res1@mpMinLatF  = -30.
    res1@mpMaxLatF  = 30.
    res1@mpMinLonF  = 120.
    res1@mpMaxLonF  = 290.
    res1@mpCenterLonF   = 205.
    ;res1@gsnSpreadColors = True
    ;res1@gsnSpreadColorStart = -60
    res1@lbOrientation       = "vertical"
    ;res1@cnLevelSelectionMode = "ExplicitLevels"
    ;res1@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)

    dum = new(24,graphic)
    dum2= new(24,graphic)
    ypts4 = (/5,5,-5,-5,5/)
    xpts4 = (/160,210,210,160,160/)
    ypts3 = (/5,5,-5,-5,5/)
    xpts3 = (/210,270,270,210,210/)

    resl4 = True
    resl4@gsLineColor = "red"
    resl4@gsLineThicknessF = 3.
;   resl4@gsLineLabelString = "Nino 4"
    resl3 = True
    resl3@gsLineColor = "black"
    resl3@gsLineThicknessF = 3.
;   resl3@gsLineLabelString = "Nino 3"
     
    
    res1@gsnLeftString  ="GPCP"
    res1@gsnRightString ="mean="+sprintf("%3.2f", precip_obs_mon_ave)+"mm/day"
    plot1(0) = gsn_csm_contour_map_ce(wks_pre,precip_obs_mon,res1)

    do i=0,3
        dum(i) = gsn_add_polyline(wks_pre,plot1(0),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(i) = gsn_add_polyline(wks_pre,plot1(0),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do

    res1@gsnLeftString  ="CESM"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(precip_cesm_mon_ave/precip_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", precip_cesm_mon_ave)+"mm/day"
    plot1(1) = gsn_csm_contour_map_ce(wks_pre,precip_cesm_mon,res1)
    do i=0,3
        dum(4+i) = gsn_add_polyline(wks_pre,plot1(1),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(4+i) = gsn_add_polyline(wks_pre,plot1(1),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do


    res1@gsnLeftString  ="FGOALS-g2"
    res1@gsnCenterString = "RMSE="+sprintf("%3.2f", abs(precip_fgoals_mon_ave/precip_obs_mon_ave-1.))
    res1@gsnRightString ="mean="+sprintf("%3.2f", precip_fgoals_mon_ave)+"mm/day"
    plot1(2) = gsn_csm_contour_map_ce(wks_pre,precip_fgoals_mon,res1)
    do i=0,3
        dum(8+i) = gsn_add_polyline(wks_pre,plot1(2),xpts4(i:i+1),ypts4(i:i+1),resl4)
        dum2(8+i) = gsn_add_polyline(wks_pre,plot1(2),xpts3(i:i+1),ypts3(i:i+1),resl3)
    end do
   
    
    pres=True
    pres@txString           ="Tropical Pacific precipitation"
    pres@gsnPanelLabelBar   = True       ; common label bar
        ;pres@gsnFrame            = False      ; don't advance frame yet
        ;pres1@lbOrientation       = "vertical"
        ;pres@lbLabelFontHeightF  = 0.015
        ;pres@pmLabelBarWidthF = 0.85
    ;pres@cnLevelSelectionMode = "ExplicitLevels"
    ;pres@cnLevels = (/0.2,0.5,1,2,3,4,5,6,7,8,9,10,12,14,17/)
    gsn_panel(wks_pre,plot1,(/3,1/),pres)
end