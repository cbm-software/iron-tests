gfx r n Example.part0.exnode
gfx r e Example.part0.exelem

gfx def faces

gfx def field DeformedGeometry composite Dependent.1 Dependent.2 Dependent.3
gfx def field Displacement add fields DeformedGeometry Undeformed scale_factors 1 -1
gfx def field DisplacementX component Displacement.1
gfx def field DisplacementY component Displacement.2
gfx def field DisplacementZ component Displacement.3
gfx define field Pressure component Dependent.4
gfx def field deludelnX component Derivative.1
gfx def field deludelnY component Derivative.2
gfx def field deludelnZ component Derivative.3

gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3
gfx modify g_element "/" general clear
#gfx modify g_element "/" node_points subgroup Region coordinate Undeformed LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label Dependent select_on material black selected_material default_selected

# create spectrum and color bar
gfx create spectrum displacement_spectrum
gfx create spectrum scalarfield_spectrum

# modify surface and point data representation
gfx modify g_element "/" surfaces coordinate Undeformed tessellation default LOCAL select_on material black data Pressure spectrum scalarfield_spectrum selected_material default_selected render_shaded
gfx modify g_element "/" surfaces coordinate DeformedGeometry tessellation default LOCAL select_on material black data DisplacementX spectrum displacement_spectrum selected_material default_selected render_shaded

gfx modify spectrum displacement_spectrum autorange
gfx modify spectrum scalarfield_spectrum autorange
gfx create colour_bar spectrum displacement_spectrum label_material black

gfx modify g_element "/" point coordinate Undeformed NORMALISED_WINDOW_FIT_LEFT glyph colour_bar general size "1*1*1" centre 0,0,0 font default select_on material black selected_material copper

# modify window and change to white background
gfx create window 1 double_buffer;
gfx modify window 1 image scene default light_model default;
gfx modify window 1 image add_light default;
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 770 height 640;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 2.52759 -0.986566 1.39318 interest_point 0.772664 0.485626 0.354958 up_vector -0.324396 0.255576 0.91074 view_angle 89.1034 near_clipping_plane 0.0251496 far_clipping_plane 8.98759 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

# color bar with black font
gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on material black selected_material black normalised_window_fit_left

# open scene editor
gfx edit scene
