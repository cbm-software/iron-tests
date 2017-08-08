#Read in the sequence of nodal positions.
gfx read node results/current_run/3D_MeshRefinementLevel_001/Example.part0.exnode time 0
for ($i = 100 ; $i < 1001 ; $i=$i+100)
{
    $filename = sprintf("results/current_run/3D_MeshRefinementLevel_001/Example_%d.part0.exnode", $i);
    gfx read node "$filename" time $i;
}

#Read in the element description
gfx read elements results/reference/iron//3D_MeshRefinementLevel_001/Example.part0.exelem;

gfx define field Velocity composite VectorField.1 VectorField.2 VectorField.3
gfx define field VectorFieldMagnitude magnitude field Velocity

gfx define field VelocityX component VectorField.1
gfx define field VelocityY component VectorField.2
gfx define field VelocityY component VectorField.3
gfx define field Pressure  component VectorField.4

gfx define faces egroup Region
gfx modify g_element Region lines coordinate Geometry select_on material default selected_material default_selected
gfx modify g_element Region node_points coordinate Geometry glyph sphere General size "0.1*0.1*0.1" centre 0,0,0 font default select_on material default selected_material default_selected

gfx cre win 1
gfx draw axes
gfx edit scene
gfx cre mat copper ambient 1 0.2 0 diffuse 0.6 0.3 0 specular 0.7 0.7 0.5 shininess 0.3
gfx modify g_element "/" general clear;

gfx create spectrum velocity_spectrum

gfx modify g_element "/" surfaces coordinate Geometry tessellation default LOCAL select_on material black spectrum velocity_spectrum selected_material default_selected render_wireframe;
gfx modify spectrum velocity_spectrum linear reverse range 0 1.0 extend_above extend_below rainbow colour_range 0 1 component 1
gfx modify g_element "/" node_points subgroup Region coordinate Geometry LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label Geometry select_on material black selected_material default_selected;
gfx modify g_element "/" surfaces coordinate Geometry tessellation default LOCAL select_on material black data VectorFieldMagnitude spectrum velocity_spectrum selected_material default_selected render_shaded;

gfx create colour_bar spectrum velocity_spectrum label_material black

gfx modify g_element "/" point glyph colour_bar general size "1*1*1" centre 0,0,0 select_on  normalised_window_fit_left;

gfx modify window 1 image view_all
gfx modify window 1 background colour 1 1 1
gfx modify window 1 layout simple ortho_axes z -y eye_spacing 0.25 width 878 height 574;
gfx modify window 1 set current_pane 1;
gfx modify window 1 background colour 1 1 1 texture none;
gfx modify window 1 view parallel eye_point 0.5 0.520833 3.06574 interest_point 0.5 0.520833 0 up_vector -0 1 0 view_angle 40 near_clipping_plane 0.0306574 far_clipping_plane 10.9559 relative_viewport ndc_placement -1 1 2 2 viewport_coordinates 0 0 1 1;
gfx modify window 1 set transform_tool current_pane 1 std_view_angle 40 normal_lines no_antialias depth_of_field 0.0 fast_transparency blend_normal;

gfx set time 0

gfx print postscript file doc/figures/current_run.eps
exit
