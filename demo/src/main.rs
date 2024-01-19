use std::fs::File;
use std::io::{prelude::*, BufReader};
use raylib::prelude::*;

#[derive(Copy, Clone, Debug)]
struct Vec3D {
    x: f32,
    y: f32,
    z: f32,
    w: f32
}

struct Triangle {
    p: [Vec3D; 3],
    color: Option<Color>
}

struct Mesh {
    tris: Vec<Triangle>
}

#[derive(Debug)]
struct Mat4x4 {
    m: [[f32; 4]; 4]
}

impl Vec3D {
    fn new(x: f32, y: f32, z: f32, w: f32) -> Self {
        Self { x, y, z, w }
    }

    fn vector_add(v1: &Self, v2: &Self) -> Self {
        Self::new(
            v1.x + v2.x,
            v1.y + v2.y,
            v1.z + v2.z,
            0.0
        )
    }

    fn vector_sub(v1: &Self, v2: &Self) -> Self {
        Self::new(
            v1.x - v2.x,
            v1.y - v2.y,
            v1.z - v2.z,
            0.0
        )
    }

    fn vector_mul(v1: &Self, k: f32) -> Self {
        Self::new(
            v1.x * k,
            v1.y * k,
            v1.z * k,
            0.0
        )
    }

    fn vector_div(v1: &Self, k: f32) -> Self {
        Self::new(
            v1.x / k,
            v1.y / k,
            v1.z / k,
            0.0
        )
    }

    fn vector_dotproduct(v1: &Self, v2: &Self) -> f32 {
        v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
    }

    fn vector_crossproduct(v1: &Self, v2: &Self) -> Self {
        let x = v1.y * v2.z - v1.z * v2.y;
        let y = v1.z * v2.x - v1.x * v2.z;
        let z = v1.x * v2.y - v1.y * v2.x;

        Self::new(x, y, z, 1.0)
    }

    fn vector_length(v: &Self) -> f32 {
        Self::vector_dotproduct(v, v).sqrt()
    }

    fn vector_normalize(v: &Self) -> Self {
        let l = Self::vector_length(v);
        Self::new(
            v.x / l,
            v.y / l,
            v.z / l,
            1.0
        )
    }

    fn matrix_multiply_vector(m: &Mat4x4, i: &Self) -> Self {
        let x = i.x * m.m[0][0] + i.y * m.m[1][0] + i.z * m.m[2][0] + i.w * m.m[3][0];
        let y = i.x * m.m[0][1] + i.y * m.m[1][1] + i.z * m.m[2][1] + i.w * m.m[3][1];
        let z = i.x * m.m[0][2] + i.y * m.m[1][2] + i.z * m.m[2][2] + i.w * m.m[3][2];
        let w = i.x * m.m[0][3] + i.y * m.m[1][3] + i.z * m.m[2][3] + i.w * m.m[3][3];

        Self::new(x, y, z, w)
    }
}

impl Mat4x4 {
    fn matrix_make_identity() -> Self {
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = 1.0;
        matrix.m[1][1] = 1.0;
        matrix.m[2][2] = 1.0;
        matrix.m[3][3] = 1.0;

        matrix
    }

    fn matrix_make_rotation_x(f_angle_rad: f32) -> Self {
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = 1.0;
        matrix.m[1][1] = f_angle_rad.cos();
        matrix.m[1][2] = f_angle_rad.sin();
        matrix.m[2][1] = -f_angle_rad.sin();
        matrix.m[2][2] = f_angle_rad.cos();
        matrix.m[3][3] = 1.0;

        matrix
    }

    fn matrix_make_rotation_y(f_angle_rad: f32) -> Self {
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = f_angle_rad.cos();
        matrix.m[0][2] = f_angle_rad.sin();
        matrix.m[2][0] = -f_angle_rad.sin();
        matrix.m[1][1] = 1.0;
        matrix.m[2][2] = f_angle_rad.cos();
        matrix.m[3][3] = 1.0;

        matrix
    }

    fn matrix_make_rotation_z(f_angle_rad: f32) -> Self {
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = f_angle_rad.cos();
        matrix.m[0][1] = f_angle_rad.sin();
        matrix.m[1][0] = -f_angle_rad.sin();
        matrix.m[1][1] = f_angle_rad.cos();
        matrix.m[2][2] = 1.0;
        matrix.m[3][3] = 1.0;

        matrix
    }

    fn matrix_make_translation(x: f32, y: f32, z: f32) -> Self {
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = 1.0;
        matrix.m[1][1] = 1.0;
        matrix.m[2][2] = 1.0;
        matrix.m[3][3] = 1.0;
        matrix.m[3][0] = x;
        matrix.m[3][1] = y;
        matrix.m[3][2] = z;

        matrix
    }

    fn matrix_make_projection(f_fov_degrees: f32, f_aspect_ratio: f32, f_near: f32, f_far: f32) -> Self {
        let f_fov_rad = 1.0 / (f_fov_degrees * 0.5 / 180.0 * 3.14159).tan();
        let mut matrix = Mat4x4::default();
        matrix.m[0][0] = f_aspect_ratio * f_fov_rad;
        matrix.m[1][1] = f_fov_rad;
        matrix.m[2][2] = f_far / (f_far - f_near);
        matrix.m[3][2] = (-f_far * f_near) / (f_far - f_near);
        matrix.m[2][3] = 1.0;
        matrix.m[3][3] = 0.0;

        matrix
    }

    fn matrix_multiply_matrix(m1: &Self, m2: &Self) -> Self {
        let mut matrix = Mat4x4::default();
        for c in 0..4 {
            for r in 0..4 {
                matrix.m[r][c] =
                    m1.m[r][0] * m2.m[0][c] +
                    m1.m[r][1] * m2.m[1][c] +
                    m1.m[r][2] * m2.m[2][c] +
                    m1.m[r][3] * m2.m[3][c];
            }
        }

        matrix
    }
}

impl Triangle {
    fn new(p: [Vec3D; 3]) -> Self {
        Self { p, color: None }
    }
}

impl Mesh {
    fn new(tris: Vec<Triangle>) -> Self {
        Self { tris }
    }

    fn load_from_object_file(&mut self, filename: &str) -> bool {
        let full_filename = format!("3d_models/{}", filename);

        let file = File::open(full_filename).unwrap();
        let reader = BufReader::new(file);

        let mut verts: Vec<Vec3D> = Vec::new();

        for line in reader.lines() {
            let un_line = line.unwrap();
            let junk = un_line.chars().nth(0).unwrap();

            match junk {
                'v' => {
                    let line_vec: Vec<&str> = un_line.split(' ').collect();

                    let x = line_vec.get(1).unwrap().to_owned().parse::<f32>().unwrap();
                    let y = line_vec.get(2).unwrap().to_owned().parse::<f32>().unwrap();
                    let z = line_vec.get(3).unwrap().to_owned().parse::<f32>().unwrap();

                    let v = Vec3D::new(x, y, z, 1.0);
                    verts.push(v);
                },
                'f' => {
                    let mut f: [i32; 3] = [0, 0, 0];
                    let line_vec: Vec<&str> = un_line.split(' ').collect();

                    f[0] = line_vec.get(1).unwrap().to_owned().parse::<i32>().unwrap();
                    f[1] = line_vec.get(2).unwrap().to_owned().parse::<i32>().unwrap();
                    f[2] = line_vec.get(3).unwrap().to_owned().parse::<i32>().unwrap();

                    let t = Triangle::new([
                        verts[(f[0] - 1) as usize],
                        verts[(f[1] - 1) as usize],
                        verts[(f[2] - 1) as usize]
                    ]);
                    self.tris.push(t);
                },
                _ => {}
            }
        }

        true
    }
}

impl Default for Mat4x4 {
    fn default() -> Self {
        Self {
            m: [[0.0; 4]; 4]
        }
    }
}

impl Default for Vec3D {
    fn default() -> Self {
        Self { x: 0.0, y: 0.0, z: 0.0, w: 1.0 }
    }
}

impl Default for Triangle {
    fn default() -> Self {
        let p0 = Vec3D::default();
        let p1 = Vec3D::default();
        let p2 = Vec3D::default();

        Self { p: [p0, p1, p2], color: None }
    }
}

impl Default for Mesh {
    fn default() -> Self {
        Self { tris: vec!() }
    }
}

fn main() {
    let mut mesh = Mesh::default();

    mesh.load_from_object_file("VideoShip.obj");

    let v_camera = Vec3D::default();
    let f_near: f32 = 0.1;
    let f_far: f32 = 1000.0;
    let f_fov: f32 = 90.0;
    let f_aspect_ratio: f32 = 640.0 / 480.0;
    let mut f_theta: f64 = 0.0;

    let mat_proj = Mat4x4::matrix_make_projection(
        f_fov,
        f_aspect_ratio,
        f_near,
        f_far
    );

    let mut mat_rot_z: Mat4x4;
    let mut mat_rot_x: Mat4x4;
    let mut mat_translation: Mat4x4;
    let mut mat_world: Mat4x4;

    let (mut rl, thread) = raylib::init()
        .size(800, 800)
        .title("3D demo")
        .build();

    while !rl.window_should_close() {
        let mut d = rl.begin_drawing(&thread);
        d.clear_background(Color::BLACK);

        f_theta = d.get_time() + 100.0;

        mat_rot_z = Mat4x4::matrix_make_rotation_z((f_theta * 0.5) as f32);
        mat_rot_x = Mat4x4::matrix_make_rotation_x(f_theta as f32);
        mat_translation = Mat4x4::matrix_make_translation(0.0, 0.0, 16.0);

        mat_world = Mat4x4::matrix_multiply_matrix(&mat_rot_x, &mat_rot_z);
        mat_world = Mat4x4::matrix_multiply_matrix(&mat_world, &mat_translation);

        let mut vec_triangles_to_raster: Vec<Triangle> = vec![];

        for tri in &mesh.tris {
            let mut tri_projected = Triangle::default();
            let mut tri_transformed = Triangle::default();

            // Draw triangles
            tri_transformed.p[0] = Vec3D::matrix_multiply_vector(&mat_world, &tri.p[0]);
            tri_transformed.p[1] = Vec3D::matrix_multiply_vector(&mat_world, &tri.p[1]);
            tri_transformed.p[2] = Vec3D::matrix_multiply_vector(&mat_world, &tri.p[2]);

            let line1 = Vec3D::vector_sub(&tri_transformed.p[1], &tri_transformed.p[0]);
            let line2 = Vec3D::vector_sub(&tri_transformed.p[2], &tri_transformed.p[0]);

            let mut normal = Vec3D::vector_crossproduct(&line1, &line2);
            normal = Vec3D::vector_normalize(&normal);

            let v_camera_ray = Vec3D::vector_sub(&tri_transformed.p[0], &v_camera);

            if Vec3D::vector_dotproduct(&normal, &v_camera_ray) < 0.0 {
                // Illumination
                let mut light_direction = Vec3D { x: 0.0, y: 0.0, z: -1.0, w: 1.0 };
                light_direction = Vec3D::vector_normalize(&light_direction);

                let dp = (0.1_f32).max(Vec3D::vector_dotproduct(&light_direction, &normal));
                let c = get_color(dp);
                tri_projected.color = Some(c);

                // Project triangles from 3D --> 2D
                tri_projected.p[0] = Vec3D::matrix_multiply_vector(&mat_proj, &tri_transformed.p[0]);
                tri_projected.p[1] = Vec3D::matrix_multiply_vector(&mat_proj, &tri_transformed.p[1]);
                tri_projected.p[2] = Vec3D::matrix_multiply_vector(&mat_proj, &tri_transformed.p[2]);

                tri_projected.p[0] = Vec3D::vector_div(&tri_projected.p[0], tri_projected.p[0].w);
                tri_projected.p[1] = Vec3D::vector_div(&tri_projected.p[1], tri_projected.p[1].w);
                tri_projected.p[2] = Vec3D::vector_div(&tri_projected.p[2], tri_projected.p[2].w);

                // Scale into view (?)
                let v_offset_view = Vec3D::new(1.0, 1.0, 0.0, 1.0);
                tri_projected.p[0] = Vec3D::vector_add(&tri_projected.p[0], &v_offset_view);
                tri_projected.p[1] = Vec3D::vector_add(&tri_projected.p[1], &v_offset_view);
                tri_projected.p[2] = Vec3D::vector_add(&tri_projected.p[2], &v_offset_view);

                tri_projected.p[0].x *= 0.5 * 800.0;
                tri_projected.p[0].y *= 0.5 * 800.0;

                tri_projected.p[1].x *= 0.5 * 800.0;
                tri_projected.p[1].y *= 0.5 * 800.0;

                tri_projected.p[2].x *= 0.5 * 800.0;
                tri_projected.p[2].y *= 0.5 * 800.0;

                // Store triangle for sorting
                vec_triangles_to_raster.push(tri_projected);
            }
        }

        // Sort triangles from back to front
        vec_triangles_to_raster.sort_by(|t1, t2| {
            let z1 = (t1.p[0].z + t1.p[1].z + t1.p[2].z) / 3.0;
            let z2 = (t2.p[0].z + t2.p[1].z + t2.p[2].z) / 3.0;

            z2.partial_cmp(&z1).unwrap()
        });

        for tri_projected in vec_triangles_to_raster {
            // Rasterize triangle
            draw_triangle(
                &mut d,
                tri_projected.p[0].x,
                tri_projected.p[0].y,
                tri_projected.p[1].x,
                tri_projected.p[1].y,
                tri_projected.p[2].x,
                tri_projected.p[2].y,
                tri_projected.color
            );
        }
    }
}

fn multiply_matrix_vector(input: &Vec3D, output: &mut Vec3D, mat_proj: &Mat4x4) {
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

fn get_color(lum: f32) -> Color {
    let pixel_bw = (13.0*lum) as i32;
    match pixel_bw {
        0 => Color::new(0, 0, 0, 255),
        1 => Color::new(19, 19, 19, 255),
        2 => Color::new(38, 38, 38, 255),
        3 => Color::new(57, 57, 57, 255),
        4 => Color::new(76, 76, 76, 255),
        5 => Color::new(95, 95, 95, 255),
        6 => Color::new(114, 114, 114, 255),
        7 => Color::new(133, 133, 133, 255),
        8 => Color::new(152, 152, 152, 255),
        9 => Color::new(171, 171, 171, 255),
        10 => Color::new(190, 190, 190, 255),
        11 => Color::new(209, 209, 209, 255),
        12 => Color::new(228, 228, 228, 255),
        _ => Color::new(0, 0, 0, 255)
    }
}

fn draw_triangle(d: &mut RaylibDrawHandle, x1: f32, y1: f32, x2: f32, y2: f32, x3: f32, y3: f32, color: Option<Color>) {
    match color {
        Some(c) => {
            let t1 = Vector2::new(x1, y1);
            let t2 = Vector2::new(x2, y2);
            let t3 = Vector2::new(x3, y3);

            d.draw_triangle(t1, t2, t3, c);

            // d.draw_line(x1 as i32, y1 as i32, x2 as i32, y2 as i32, Color::BLACK);
            // d.draw_line(x2 as i32, y2 as i32, x3 as i32, y3 as i32, Color::BLACK);
            // d.draw_line(x3 as i32, y3 as i32, x1 as i32, y1 as i32, Color::BLACK);
        },
        None => {
            eprintln!("Triangle doesn't have a color value.");
            std::process::exit(1);
        }
    }
}
