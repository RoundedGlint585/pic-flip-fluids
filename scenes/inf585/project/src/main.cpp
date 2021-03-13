/**
	Introduction to use the VCL library
*/


#include "vcl/vcl.hpp"
#include <iostream>
#include "Particles.h"

// Add vcl namespace within the current one - Allows to use function from vcl library without explicitely preceeding their name with vcl::
using namespace vcl;


// ****************************************** //
// Structures associated to the current scene
//   In general these structures will be pre-coded 
//   for you in the exercises
// ****************************************** //

// Structure storing the variables used in the GUI interface
struct gui_parameters {
	bool display_color = true;
	bool display_particles = true;
	bool display_radius = false;
};


// Structure storing user-related interaction data and GUI parameter
struct user_interaction_parameters {
	vec2 mouse_prev;
	timer_fps fps_record;
	gui_parameters gui;
	bool cursor_on_gui;
};
user_interaction_parameters user;

// Structure storing the global variable of the 3D scene
//   can be use to send uniform parameter when displaying a shape
struct scene_environment
{
	camera_around_center camera; // A camera looking at, and rotating around, a specific center position
	mat4 projection;             // The perspective projection matrix
	vec3 light;                  // Position of the light in the scene
};
scene_environment scene;  // (declaration of scene as a global variable)


// ****************************************** //
// Functions signatures
// ****************************************** //

// Callback functions
//   Functions called when a corresponding event is received by GLFW (mouse move, keyboard, etc).
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos);
void mouse_click_callback(GLFWwindow* window, int button, int action, int mods);
void window_size_callback(GLFWwindow* window, int width, int height);

// Specific functions:
//    Initialize the data of this scene - executed once before the start of the animation loop.
void initialize_data();   
void display_interface();
//    Display calls - executed a each frame
void display_scene(float current_time);
//    Display the GUI widgets
void update_field_color(MACGrid& grid, Particles& particles);



// ****************************************** //
// Declaration of Global variables
// ****************************************** //
timer_basic timer;
MACGrid grid(25, 25, 0.04);
mesh_drawable field_quad;
Particles particles(4, grid);
mesh_drawable sphere_particle;
curve_drawable curve_visual;
mesh_drawable box;


// ****************************************** //
// Functions definitions
// ****************************************** //


// Main function with creation of the scene and animation loop
int main(int, char* argv[])
{
	std::cout << "Run " << argv[0] << std::endl;

	// create GLFW window and initialize OpenGL
	GLFWwindow* window = create_window(1280,1024); 
	window_size_callback(window, 1280, 1024);
	std::cout << opengl_info_display() << std::endl;;

	imgui_init(window); // Initialize GUI library

	// Set GLFW callback functions
	glfwSetCursorPosCallback(window, mouse_move_callback); 
	glfwSetWindowSizeCallback(window, window_size_callback);
	glfwSetMouseButtonCallback(window, mouse_click_callback);
	
	std::cout<<"Initialize data ..."<<std::endl;
	initialize_data();

	std::cout<<"Start animation loop ..."<<std::endl;
	user.fps_record.start();
	timer.start();
	glEnable(GL_DEPTH_TEST);
	while (!glfwWindowShouldClose(window))
	{
		scene.light = scene.camera.position();
		user.fps_record.update();
		timer.update(); // update the time at this current frame

		// Clear screen
		glClearColor(1.0f, 1.0f, 1.0f, 1.0f);
		glClear(GL_COLOR_BUFFER_BIT);
		glClear(GL_DEPTH_BUFFER_BIT);

		// Create GUI interface for the current frame
		imgui_create_frame();
		if (user.fps_record.event) {
			std::string const title = "VCL Display - " + str(user.fps_record.fps) + " fps";
			glfwSetWindowTitle(window, title.c_str());
		}
		ImGui::Begin("GUI",NULL,ImGuiWindowFlags_AlwaysAutoResize);
		user.cursor_on_gui = ImGui::IsAnyWindowFocused();

		float dt = 0.02f * timer.scale;
        //dt = grid.CFL()/3;
		particles.step(dt);
		
		// Set the GUI interface (widgets: buttons, checkbox, sliders, etc)
		display_interface();

		// Display the objects of the scene
		display_scene(timer.t);


		// Display GUI
		ImGui::End();
		imgui_render_frame(window);

		// Swap buffer and handle GLFW events
		glfwSwapBuffers(window);
		glfwPollEvents();
	}


	imgui_cleanup();
	glfwDestroyWindow(window);
	glfwTerminate();

	return 0;
}

void initialize_pic() {
	//
}

void initialize_data()
{
	// Load and set the common shaders
	// *************************************** //

	//   - Shader used to display meshes
	GLuint const shader_mesh = opengl_create_shader_program(opengl_shader_preset("mesh_vertex"), opengl_shader_preset("mesh_fragment"));
	//   - Shader used to display constant color (ex. for curves)
	GLuint const shader_uniform_color = opengl_create_shader_program(opengl_shader_preset("single_color_vertex"), opengl_shader_preset("single_color_fragment"));
	//   - Default white texture
	GLuint const texture_white = opengl_texture_to_gpu(image_raw{1,1,image_color_type::rgba,{255,255,255,255}});

	// Set default shader and texture to drawable mesh
	mesh_drawable::default_shader = shader_mesh;
	mesh_drawable::default_texture = texture_white;
	curve_drawable::default_shader = shader_uniform_color;
	segments_drawable::default_shader = shader_uniform_color;

	// Set the initial position of the camera
	// *************************************** //

	vec3 const camera_position = { 0,0,3.0f };        // position of the camera in space
	vec3 const camera_target_position = {0,0,0}; // position the camera is looking at / point around which the camera rotates
	vec3 const up = {0,1,0};                     // approximated "up" vector of the camera
	scene.camera.look_at(camera_position, camera_target_position, up);

	// Prepare the objects visible in the scene
	// *************************************** //
	field_quad = mesh_drawable(mesh_primitive_quadrangle({ -1,-1,0 }, { 1,-1,0 }, { 1,1,0 }, { -1,1,0 }));
	field_quad.shading.phong = { 1,0,0 };
	// field_quad.texture = opengl_texture_to_gpu(grid.getDensity);

	initialize_pic();
	sphere_particle = mesh_drawable(mesh_primitive_sphere());
	sphere_particle.transform.scale = .01f;
	box = mesh_drawable(mesh_primitive_grid());
	curve_visual.color = { 1,0,0 };
	curve_visual = curve_drawable(curve_primitive_circle());
}



void display_scene(float time)
{
	auto positions = particles.getParticlePositions();
	if (user.gui.display_particles) {
		for (auto & position : positions) {
			vec3 const& p = vec3 (position, 0);
			sphere_particle.transform.translate = p;
			draw(sphere_particle, scene);
			draw(box, scene);
		}
	}

	if (user.gui.display_radius) {
		curve_visual.transform.scale = 0.12f;  // Influence distance of a particle (size of the kernel), I took the same as SPH for now
		for (size_t k = 0; k < particles.getParticlePositions().size(); k++) {
			curve_visual.transform.translate = vec3(positions[k], 0);
			draw(curve_visual, scene);
		}
	}

	if (user.gui.display_color) {
		update_field_color(grid, particles);
		// opengl_update_texture_gpu(field_quad.texture, grid);
		draw(field_quad, scene);
	}
}

// Display the GUI
void display_interface()
{
	ImGui::SliderFloat("Timer scale", &timer.scale, 0.01f, 4.0f, "%0.2f");

	bool const restart = ImGui::Button("Restart");
	if (restart)
		initialize_pic();

	ImGui::Checkbox("Color", &user.gui.display_color);
	ImGui::Checkbox("Particles", &user.gui.display_particles);
	ImGui::Checkbox("Radius", &user.gui.display_radius);

}

// Function called every time the screen is resized
void window_size_callback(GLFWwindow* , int width, int height)
{	
	glViewport(0, 0, width, height); // The image is displayed on the entire window
	float const aspect = width / static_cast<float>(height); // Aspect ratio of the window

	// Generate the perspective matrix for this aspect ratio
	float const field_of_view = 50.0f*pi/180.0f; // the angle of the field of view
	float const z_near = 0.1f;  // closest visible point
	float const z_far = 100.0f; // furthest visible point
	scene.projection = projection_perspective(field_of_view, aspect, z_near, z_far); 
}

// Function called every time the mouse is moved
void mouse_move_callback(GLFWwindow* window, double xpos, double ypos)
{
	vec2 const  p1 = glfw_get_mouse_cursor(window, xpos, ypos);
	vec2 const& p0 = user.mouse_prev;

	glfw_state state = glfw_current_state(window);
	

	// Handle camera rotation
	auto& camera = scene.camera;
	if(!user.cursor_on_gui){
		if(state.mouse_click_left && !state.key_ctrl)
			scene.camera.manipulator_rotate_trackball(p0, p1);
		if(state.mouse_click_left && state.key_ctrl)
			camera.manipulator_translate_in_plane(p1-p0);
		if(state.mouse_click_right)
			camera.manipulator_scale_distance_to_center( (p1-p0).y );
	}

	user.mouse_prev = p1;
}

void mouse_click_callback(GLFWwindow* window, int button, int action, int mods)
{
	ImGui::SetWindowFocus(nullptr);
}

// Uniform data used when displaying an object in this scene
void opengl_uniform(GLuint shader, scene_environment const& current_scene)
{
	opengl_uniform(shader, "projection", current_scene.projection);
	opengl_uniform(shader, "view", scene.camera.matrix_view());
	opengl_uniform(shader, "light", scene.light, false);
}

void update_field_color(MACGrid& grid, Particles& particles)
{
	// field.fill({ 1,1,1 });
	float const d = 0.1f;
	int Nf = grid.xCellCount;
	for (int kx = 0; kx < Nf; ++kx) {
		for (int ky = 0; ky < grid.yCellCount; ++ky) {
			float f = 0.0f;
			vec3 const p0 = { 2.0f * (kx / (Nf - 1.0f) - 0.5f), 2.0f * (ky / (Nf - 1.0f) - 0.5f), 0.0f };
			auto positions = particles.getParticlePositions();
			for (size_t k = 0; k < particles.getParticlePositions().size(); ++k) {
				vec3 const& pi = vec3(positions[k], 0.0f);
				float const r = norm(pi - p0) / d;
				f += 0.25f * std::exp(-r * r);
			}
			// field(kx, Nf - 1 - ky) = vec3(clamp(1 - f, 0, 1), clamp(1 - f, 0, 1), 1);
		}
	}
}

