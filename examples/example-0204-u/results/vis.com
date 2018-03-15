gfx r n Example.part0.exnode
gfx r e Example.part0.exelem

gfx def faces

gfx def field DeformedGeometry composite Dependent.1 Dependent.2 Dependent.3
gfx def field Displacement add fields DeformedGeometry Undeformed scale_factors 1 -1
gfx def field DisplacementX component Displacement.1
gfx def field DisplacementY component Displacement.2
gfx def field DisplacementZ component Displacement.3
gfx def field DisplacementM magnitude field Displacement
gfx def field Pressure component Dependent.4
gfx def field deludelnX component Derivative.1
gfx def field deludelnY component Derivative.2
gfx def field deludelnZ component Derivative.3
gfx def field PositionX component Dependent.1
gfx def field PositionY component Dependent.2
gfx def field PositionZ component Dependent.3

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3
gfx modify g_element "/" general clear

# create spectrum
gfx create spectrum displacementx_spectrum
gfx create spectrum displacementy_spectrum
gfx create spectrum displacementz_spectrum
gfx create spectrum displacementm_spectrum
gfx create spectrum pressure_spectrum

# modify surface and point data representation
gfx modify g_element "/" surfaces coordinate Undeformed tessellation default LOCAL select_on material tissue selected_material default_selected render_shaded

# modify spectrum
gfx modify spectrum displacementx_spectrum autorange
gfx modify spectrum displacementy_spectrum autorange
gfx modify spectrum displacementz_spectrum autorange
gfx modify spectrum displacementm_spectrum autorange

# modify window and change to white background
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 459 height 640;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 584.852 -257.246 -1227.06 interest_point 105.153 145.168 -1510.86 up_vector -0.324398 0.255578 0.910739 view_angle 26.6957 near_clipping_plane 6.8745 far_clipping_plane 2456.71 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

# open scene editor
gfx edit scene

#gfx print postscript file ../../../doc/figures/undeformed_geometry.eps
gfx print postscript file ../../../doc/figures/current_run-undeformed_geometry.eps

############## plot displacement components, magnitude and pressure
# first remove undeformed domain
gfx modify g_element "/" surfaces coordinate Undeformed invisible

# pressure on deformed domain
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material tissue data Pressure spectrum pressure_spectrum selected_material default_selected render_shaded
gfx modify spectrum pressure_spectrum autorange
gfx create colour_bar spectrum pressure_spectrum label_material black
gfx modify g_element "/" point coordinate DeformedGeometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper
#gfx print postscript file ../../../doc/figures/deformed_geometry-pressure.eps
gfx print postscript file ../../../doc/figures/current_run-deformed_geometry-pressure.eps

# DisplacementM on deformed domain
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material tissue data DisplacementM spectrum displacementm_spectrum selected_material default_selected render_shaded
gfx modify spectrum displacementm_spectrum autorange
gfx create colour_bar spectrum displacementm_spectrum label_material black
gfx modify g_element "/" point coordinate DeformedGeometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper
#gfx print postscript file ../../../doc/figures/deformed_geometry-displacementM.eps
gfx print postscript file ../../../doc/figures/current_run-deformed_geometry-displacementM.eps

# DisplacementX on deformed domain
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material tissue data DisplacementX spectrum displacementx_spectrum selected_material default_selected render_shaded
gfx modify spectrum displacementx_spectrum autorange
gfx create colour_bar spectrum displacementx_spectrum label_material black
gfx modify g_element "/" point coordinate DeformedGeometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper
#gfx print postscript file ../../../doc/figures/deformed_geometry-displacementX.eps
gfx print postscript file ../../../doc/figures/current_run-deformed_geometry-displacementX.eps

# DisplacementY on deformed domain
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material tissue data DisplacementY spectrum displacementy_spectrum selected_material default_selected render_shaded
gfx modify spectrum displacementy_spectrum autorange
gfx create colour_bar spectrum displacementy_spectrum label_material black
gfx modify g_element "/" point coordinate DeformedGeometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper
#gfx print postscript file ../../../doc/figures/deformed_geometry-displacementY.eps
gfx print postscript file ../../../doc/figures/current_run-deformed_geometry-displacementY.eps

# DisplacementZ on deformed domain
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material tissue data DisplacementZ spectrum displacementz_spectrum selected_material default_selected render_shaded
gfx modify spectrum displacementz_spectrum autorange
gfx create colour_bar spectrum displacementz_spectrum label_material black
gfx modify g_element "/" point coordinate DeformedGeometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper
#gfx print postscript file ../../../doc/figures/deformed_geometry-displacementZ.eps
gfx print postscript file ../../../doc/figures/current_run-deformed_geometry-displacementZ.eps

