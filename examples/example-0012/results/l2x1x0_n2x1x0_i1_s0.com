# import node and element data
gfx read nodes results/current_run/l2x1x0_n2x1x0_i1_s0/Example.part0.exnode
gfx read elements results/current_run/l2x1x0_n2x1x0_i1_s0//Example.part0.exelem

# compute element faces
gfx define faces

# create window
gfx cre window 1 double_buffer

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3
gfx modify g_element "/" general clear
gfx modify g_element "/" node_points subgroup Region coordinate Geometry LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label ScalarField select_on material black selected_material default_selected

# create spectrum and color bar
gfx create spectrum scalarfield_spectrum
gfx create colour_bar spectrum scalarfield_spectrum label_material black

# modify surface and point data representation
gfx modify g_element "/" surfaces coordinate Geometry tessellation default LOCAL select_on material black data ScalarField spectrum scalarfield_spectrum selected_material default_selected render_shaded
gfx modify g_element "/" point coordinate Geometry NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper

# modify window and change to white background
#gfx modify window 1 image view_all
#gfx modify window 1 background colour 1 1 1

gfx modify window 1 image scene default light_model default
gfx modify window 1 image add_light default
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 770 height 640
gfx modify window 1 set current_pane 1
gfx modify window 1 background colour 1 1 1 texture none
gfx modify window 1 view parallel eye_point 4.48199 -2.61522 2.28368 interest_point 0.881887 0.332743 0.46274 up_vector -0.293807 0.216311 0.931068 view_angle 40 near_clipping_plane 0.0499671 far_clipping_plane 17.8565 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal

# color bar with black font
gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material black selected_material black normalised_window_fit_left

# open scene editor
gfx edit scene

gfx print postscript file doc/figures/current_run_l2x1x0_n2x1x0_i1_s0.eps
exit
