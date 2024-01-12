use raylib::prelude::*;

struct Vec3D {
    x: f64,
    y: f64,
    z: f64
}

struct Triangle {
    p: [Vec3D; 3]
}

struct Mesh {
    tris: Vec<Triangle>
}

impl Vec3D {
    fn new(x: f64, y: f64, z: f64) -> Self {
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
}
