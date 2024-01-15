use raylib::prelude::*;

#[derive(Copy, Clone)]
struct Vec3D {
    x: f32,
    y: f32,
    z: f32
}

struct Triangle {
    p: [Vec3D; 3]
}

struct Mesh {
    tris: Vec<Triangle>
}

struct MatProj {
    m: [[f32; 4]; 4]
}

impl Vec3D {
    fn new(x: f32, y: f32, z: f32) -> Self {
        Self { x, y, z }
    }
}

impl Triangle {
    fn new(p: [Vec3D; 3]) -> Self {
        Self { p }
    }
}

impl Mesh {
    fn new(tris: Vec<Triangle>) -> Self {
        Self { tris }
    }
}

impl Default for MatProj {
    fn default() -> Self {
        Self {
            m: [[0.0; 4]; 4]
        }
    }
}

impl Default for Vec3D {
    fn default() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0 }
    }
}

impl Default for Triangle {
    fn default() -> Self {
        let p0 = Vec3D::default();
        let p1 = Vec3D::default();
        let p2 = Vec3D::default();

        Self { p: [p0, p1, p2] }
    }
}

fn main() {
    // SOUTH
    let south_v1 = Vec3D::new(0.0, 0.0, 0.0);
    let south_v2 = Vec3D::new(0.0, 1.0, 0.0);
    let south_v3 = Vec3D::new(1.0, 1.0, 0.0);
    let south_v4 = Vec3D::new(0.0, 0.0, 0.0);
    let south_v5 = Vec3D::new(1.0, 1.0, 0.0);
    let south_v6 = Vec3D::new(1.0, 0.0, 0.0);

    let south_t1 = Triangle::new([south_v1, south_v2, south_v3]);
    let south_t2 = Triangle::new([south_v4, south_v5, south_v6]);

    // EAST
    let east_v1 = Vec3D::new(1.0, 0.0, 0.0);
    let east_v2 = Vec3D::new(1.0, 1.0, 0.0);
    let east_v3 = Vec3D::new(1.0, 1.0, 1.0);
    let east_v4 = Vec3D::new(1.0, 0.0, 0.0);
    let east_v5 = Vec3D::new(1.0, 1.0, 1.0);
    let east_v6 = Vec3D::new(1.0, 0.0, 1.0);

    let east_t1 = Triangle::new([east_v1, east_v2, east_v3]);
    let east_t2 = Triangle::new([east_v4, east_v5, east_v6]);

    // NORTH
    let north_v1 = Vec3D::new(1.0, 0.0, 1.0);
    let north_v2 = Vec3D::new(1.0, 1.0, 1.0);
    let north_v3 = Vec3D::new(0.0, 1.0, 1.0);
    let north_v4 = Vec3D::new(1.0, 0.0, 1.0);
    let north_v5 = Vec3D::new(0.0, 1.0, 1.0);
    let north_v6 = Vec3D::new(0.0, 0.0, 1.0);

    let north_t1 = Triangle::new([north_v1, north_v2, north_v3]);
    let north_t2 = Triangle::new([north_v4, north_v5, north_v6]);

    // WEST
    let west_v1 = Vec3D::new(0.0, 0.0, 1.0);
    let west_v2 = Vec3D::new(0.0, 1.0, 1.0);
    let west_v3 = Vec3D::new(0.0, 1.0, 0.0);
    let west_v4 = Vec3D::new(0.0, 0.0, 1.0);
    let west_v5 = Vec3D::new(0.0, 1.0, 0.0);
    let west_v6 = Vec3D::new(0.0, 0.0, 0.0);

    let west_t1 = Triangle::new([west_v1, west_v2, west_v3]);
    let west_t2 = Triangle::new([west_v4, west_v5, west_v6]);

    // TOP
    let top_v1 = Vec3D::new(0.0, 1.0, 0.0);
    let top_v2 = Vec3D::new(0.0, 1.0, 1.0);
    let top_v3 = Vec3D::new(1.0, 1.0, 1.0);
    let top_v4 = Vec3D::new(0.0, 1.0, 0.0);
    let top_v5 = Vec3D::new(1.0, 1.0, 1.0);
    let top_v6 = Vec3D::new(1.0, 1.0, 0.0);

    let top_t1 = Triangle::new([top_v1, top_v2, top_v3]);
    let top_t2 = Triangle::new([top_v4, top_v5, top_v6]);

    // BOTTOM
    let bottom_v1 = Vec3D::new(1.0, 0.0, 1.0);
    let bottom_v2 = Vec3D::new(0.0, 0.0, 1.0);
    let bottom_v3 = Vec3D::new(0.0, 0.0, 0.0);
    let bottom_v4 = Vec3D::new(1.0, 0.0, 1.0);
    let bottom_v5 = Vec3D::new(0.0, 0.0, 0.0);
    let bottom_v6 = Vec3D::new(1.0, 0.0, 0.0);

    let bottom_t1 = Triangle::new([bottom_v1, bottom_v2, bottom_v3]);
    let bottom_t2 = Triangle::new([bottom_v4, bottom_v5, bottom_v6]);

    let triangle_vec = vec![
        south_t1,
        south_t2,
        east_t1,
        east_t2,
        north_t1,
        north_t2,
        west_t1,
        west_t2,
        top_t1,
        top_t2,
        bottom_t1,
        bottom_t2
    ];

    let mesh = Mesh::new(triangle_vec);
    let mut mat_proj: MatProj = Default::default();

    let v_camera = Vec3D::default();
    let f_near: f32 = 0.1;
    let f_far: f32 = 1000.0;
    let f_fov: f32 = 90.0;
    let f_aspect_ratio: f32 = 640.0 / 480.0;
    let f_fov_rad: f32 = 1.0 / (f_fov * 0.5 / 180.0 * 3.14159).tan();
    let mut f_theta: f64 = 0.0;

    mat_proj.m[0][0] = f_aspect_ratio * f_fov_rad;
    mat_proj.m[1][1] = f_fov_rad;
    mat_proj.m[2][2] = f_far / (f_far - f_near);
    mat_proj.m[3][2] = (-f_far * f_near) / (f_far - f_near);
    mat_proj.m[2][3] = 1.0;
    mat_proj.m[3][3] = 0.0;

    let mut mat_rot_z: MatProj = Default::default();
    let mut mat_rot_x: MatProj = Default::default();

    let (mut rl, thread) = raylib::init()
        .size(800, 800)
        .title("3D demo")
        .build();

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);

        f_theta = d.get_time() + 100.0;

        // Rotation Z
        mat_rot_z.m[0][0] = f_theta.cos() as f32;
        mat_rot_z.m[0][1] = f_theta.sin() as f32;
        mat_rot_z.m[1][0] = -f_theta.sin() as f32;
        mat_rot_z.m[1][1] = f_theta.cos() as f32;
        mat_rot_z.m[2][2] = 1.0;
        mat_rot_z.m[3][3] = 1.0;

        // Rotation X
        mat_rot_x.m[0][0] = 1.0;
        mat_rot_x.m[1][1] = (f_theta * 0.5).cos() as f32;
        mat_rot_x.m[1][2] = (f_theta * 0.5).sin() as f32;
        mat_rot_x.m[2][1] = -(f_theta * 0.5).sin() as f32;
        mat_rot_x.m[2][2] = (f_theta * 0.5).cos() as f32;
        mat_rot_x.m[3][3] = 1.0;

        for tri in &mesh.tris {
            // Draw triangles
            let mut tri_projected = Triangle::default();
            let mut tri_rotated_z = Triangle::default();
            let mut tri_rotated_zx = Triangle::default();

            multiply_matrix_vector(&tri.p[0], &mut tri_rotated_z.p[0], &mat_rot_z);
            multiply_matrix_vector(&tri.p[1], &mut tri_rotated_z.p[1], &mat_rot_z);
            multiply_matrix_vector(&tri.p[2], &mut tri_rotated_z.p[2], &mat_rot_z);

            // Rotate in X-Axis
            multiply_matrix_vector(&tri_rotated_z.p[0], &mut tri_rotated_zx.p[0], &mat_rot_x);
            multiply_matrix_vector(&tri_rotated_z.p[1], &mut tri_rotated_zx.p[1], &mat_rot_x);
            multiply_matrix_vector(&tri_rotated_z.p[2], &mut tri_rotated_zx.p[2], &mat_rot_x);

            let mut tri_translated = Triangle::new([
                tri_rotated_zx.p[0],
                tri_rotated_zx.p[1],
                tri_rotated_zx.p[2]
            ]);

            // Offset into the screen
            tri_translated.p[0].z = tri_rotated_zx.p[0].z + 3.0;
            tri_translated.p[1].z = tri_rotated_zx.p[1].z + 3.0;
            tri_translated.p[2].z = tri_rotated_zx.p[2].z + 3.0;

            let mut normal = Vec3D::default();
            let mut line1 = Vec3D::default();
            let mut line2 = Vec3D::default();

            line1.x = tri_translated.p[1].x - tri_translated.p[0].x;
            line1.y = tri_translated.p[1].y - tri_translated.p[0].y;
            line1.z = tri_translated.p[1].z - tri_translated.p[0].z;

            line2.x = tri_translated.p[2].x - tri_translated.p[0].x;
            line2.y = tri_translated.p[2].y - tri_translated.p[0].y;
            line2.z = tri_translated.p[2].z - tri_translated.p[0].z;

            normal.x = line1.y * line2.z - line1.z * line2.y;
            normal.y = line1.z * line2.x - line1.x * line2.z;
            normal.z = line1.x * line2.y - line1.y * line2.x;

            let l = (normal.x*normal.x + normal.y*normal.y + normal.z*normal.z).sqrt();
            normal.x /= l; normal.y /= l; normal.z /= l;

            // if normal.z < 0.0 {
            if (
                normal.x * (tri_translated.p[0].x - v_camera.x) +
                normal.y * (tri_translated.p[0].y - v_camera.y) +
                normal.z * (tri_translated.p[0].z - v_camera.z)
            ) < 0.0 {
                // Project triangles from 3D --> 2D
                multiply_matrix_vector(&tri_translated.p[0], &mut tri_projected.p[0], &mat_proj);
                multiply_matrix_vector(&tri_translated.p[1], &mut tri_projected.p[1], &mat_proj);
                multiply_matrix_vector(&tri_translated.p[2], &mut tri_projected.p[2], &mat_proj);

                // Scale into view (?)
                tri_projected.p[0].x += 1.0; tri_projected.p[0].y += 1.0;
                tri_projected.p[1].x += 1.0; tri_projected.p[1].y += 1.0;
                tri_projected.p[2].x += 1.0; tri_projected.p[2].y += 1.0;

                tri_projected.p[0].x *= 0.5 * 800.0;
                tri_projected.p[0].y *= 0.5 * 800.0;

                tri_projected.p[1].x *= 0.5 * 800.0;
                tri_projected.p[1].y *= 0.5 * 800.0;

                tri_projected.p[2].x *= 0.5 * 800.0;
                tri_projected.p[2].y *= 0.5 * 800.0;

                draw_triangle(
                    &mut d,
                    tri_projected.p[0].x,
                    tri_projected.p[0].y,
                    tri_projected.p[1].x,
                    tri_projected.p[1].y,
                    tri_projected.p[2].x,
                    tri_projected.p[2].y,
                    Color::RED
                );
            }
        }
    }
}

fn multiply_matrix_vector(input: &Vec3D, output: &mut Vec3D, mat_proj: &MatProj) {
    output.x = input.x * mat_proj.m[0][0] + input.y * mat_proj.m[1][0] + input.z * mat_proj.m[2][0] + mat_proj.m[3][0];
    output.y = input.x * mat_proj.m[0][1] + input.y * mat_proj.m[1][1] + input.z * mat_proj.m[2][1] + mat_proj.m[3][1];
    output.z = input.x * mat_proj.m[0][2] + input.y * mat_proj.m[1][2] + input.z * mat_proj.m[2][2] + mat_proj.m[3][2];
    let w: f32 = input.x * mat_proj.m[0][3] + input.y * mat_proj.m[1][3] + input.z * mat_proj.m[2][3] + mat_proj.m[3][3];

    if w != 0.0 {
        output.x /= w;
        output.y /= w;
        output.z /= w;
    }
}

fn draw_triangle(d: &mut RaylibDrawHandle, x1: f32, y1: f32, x2: f32, y2: f32, x3: f32, y3: f32, color: Color) {
    d.draw_line(x1 as i32, y1 as i32, x2 as i32, y2 as i32, color);
    d.draw_line(x2 as i32, y2 as i32, x3 as i32, y3 as i32, color);
    d.draw_line(x3 as i32, y3 as i32, x1 as i32, y1 as i32, color);
}
