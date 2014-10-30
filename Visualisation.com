gfx read node region deformed
gfx read elem region deformed

gfx read node region undef
gfx read elem region undef

gfx define faces egroup deformed
gfx define faces egroup undef

gfx modify g_element "/" general clear;
gfx modify g_element "/" point LOCAL glyph axes_solid_xyz general size "0.5*0.5*0.5" centre 0,0,0 font default select_on material default selected_material default_selected;
gfx modify g_element /deformed/ general clear;
gfx modify g_element /deformed/ lines coordinate DeformedGeometry tessellation default LOCAL select_on material default selected_material default_selected;
gfx modify g_element /deformed/ node_points coordinate DeformedGeometry LOCAL glyph point general size "1*1*1" centre 0,0,0 font default label DeformedGeometry select_on material default selected_material default_selected;
gfx modify g_element /deformed/ node_points coordinate DeformedGeometry LOCAL glyph point general size "1*1*1" centre -0.05,0,0 font default label cmiss_number select_on material red selected_material default_selected;
gfx modify g_element /undef/ general clear;
gfx modify g_element /undef/ lines coordinate Geometry tessellation default LOCAL select_on material green selected_material default_selected;

gfx cre win
