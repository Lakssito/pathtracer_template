#define _CRT_SECURE_NO_WARNINGS 1
#include <vector>
#include <cmath>
#include <random>
#include <omp.h>


#include <map>
#include <string>
#include <fstream>








#define STB_IMAGE_WRITE_IMPLEMENTATION
#include "stb_image_write.h"








#define STB_IMAGE_IMPLEMENTATION
#include "stb_image.h"








#ifndef M_PI
#define M_PI 3.14159265358979323856
#endif








static std::default_random_engine engine[32];
static std::uniform_real_distribution<double> uniform(0, 1);








double sqr(double x) { return x * x; };








class Vector {
public:
 explicit Vector(double x = 0, double y = 0, double z = 0) {
     data[0] = x;
     data[1] = y;
     data[2] = z;
 }
 double norm2() const {
     return data[0] * data[0] + data[1] * data[1] + data[2] * data[2];
 }
 double norm() const {
     return sqrt(norm2());
 }
 void normalize() {
     double n = norm();
     data[0] /= n;
     data[1] /= n;
     data[2] /= n;
 }
 double operator[](int i) const { return data[i]; };
 double& operator[](int i) { return data[i]; };
 double data[3];
};








Vector operator+(const Vector& a, const Vector& b) {
 return Vector(a[0] + b[0], a[1] + b[1], a[2] + b[2]);
}
Vector operator-(const Vector& a, const Vector& b) {
 return Vector(a[0] - b[0], a[1] - b[1], a[2] - b[2]);
}
Vector operator*(const double a, const Vector& b) {
 return Vector(a*b[0], a*b[1], a*b[2]);
}
Vector operator*(const Vector& a, const double b) {
 return Vector(a[0]*b, a[1]*b, a[2]*b);
}
Vector operator*(const Vector& a, const Vector& b) {
 return Vector(a[0]*b[0], a[1]*b[1], a[2]*b[2]);
}
Vector operator/(const Vector& a, const double b) {
 return Vector(a[0] / b, a[1] / b, a[2] / b);
}
double dot(const Vector& a, const Vector& b) {
 return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}
Vector cross(const Vector& a, const Vector& b) {
 return Vector(a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]);
}








class Ray {
public:
 Ray(const Vector& origin, const Vector& unit_direction) : O(origin), u(unit_direction) {};
 Vector O, u;
};








class Object {
public:
 Object(const Vector& albedo, bool mirror = false, bool transparent = false) : albedo(albedo), mirror(mirror), transparent(transparent) {};








 virtual bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const = 0;








 Vector albedo;
 bool mirror, transparent;
};








class Sphere : public Object {
public:
 Sphere(const Vector& center, double radius, const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent), C(center), R(radius) {};








 // returns true iif there is an intersection between the ray and the sphere
 // if there is an intersection, also computes the point of intersection P,
 // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
 // and the unit normal N
 bool intersect(const Ray& ray, Vector& P, double &t, Vector& N) const {
      // TODO (lab 1) : compute the intersection (just true/false at the begining of lab 1, then P, t and N as well)
     double b = dot(ray.u, ray.O - C);
     double delta = sqr(b) - (dot(ray.O - C, ray.O - C) - R*R);
     if (delta < 0) return false;




     // ROOTS of the quadriatci function
     double t2 = -b - sqrt(delta);
     double t1 = -b + sqrt(delta);




     // Note : t1 and t2 always positive but check for sanity




     // t2 first since t1 its always larger than t2
     if (t2 > 0) {
         t = t2;
     } else if (t1 > 0) {
         t = t1;
     } else {
         return false;
     }
     P = ray.O + t * ray.u;
     N = P - C;
     N.normalize();
     return true;
 }








 double R;
 Vector C;
};










class TriangleIndices {
public:
   TriangleIndices(int vtxi = -1, int vtxj = -1, int vtxk = -1, int ni = -1, int nj = -1, int nk = -1, int uvi = -1, int uvj = -1, int uvk = -1, int group = -1) {
       vtx[0] = vtxi; vtx[1] = vtxj; vtx[2] = vtxk;
       uv[0] = uvi; uv[1] = uvj; uv[2] = uvk;
       n[0] = ni; n[1] = nj; n[2] = nk;
       this->group = group;
   };
   int vtx[3]; // indices within the vertex coordinates array
   int uv[3];  // indices within the uv coordinates array
   int n[3];   // indices within the normals array
   int group;  // face group
};


// Class only used in labs 3 and 4






// I will provide you with an obj mesh loader (labs 3 and 4)
class TriangleMesh : public Object {
public:
   TriangleMesh(const Vector& albedo, bool mirror = false, bool transparent = false) : ::Object(albedo, mirror, transparent) {};


   // first scale and then translate the current object
   void scale_translate(double s, const Vector& t) {
       for (int i = 0; i < vertices.size(); i++) {
           vertices[i] = vertices[i] * s + t;
       }
   }


   // read an .obj file
   void readOBJ(const char* obj) {
       std::ifstream f(obj);
       if (!f) return;


       std::map<std::string, int> mtls;
       int curGroup = -1, maxGroup = -1;


       // OBJ indices are 1-based and can be negative (relative), this normalizes them
       auto resolveIdx = [](int i, int size) {
           return i < 0 ? size + i : i - 1;
       };


       auto setFaceVerts = [&](TriangleIndices& t, int i0, int i1, int i2) {
           t.vtx[0] = resolveIdx(i0, vertices.size());
           t.vtx[1] = resolveIdx(i1, vertices.size());
           t.vtx[2] = resolveIdx(i2, vertices.size());
       };
       auto setFaceUVs = [&](TriangleIndices& t, int j0, int j1, int j2) {
           t.uv[0] = resolveIdx(j0, uvs.size());
           t.uv[1] = resolveIdx(j1, uvs.size());
           t.uv[2] = resolveIdx(j2, uvs.size());
       };
       auto setFaceNormals = [&](TriangleIndices& t, int k0, int k1, int k2) {
           t.n[0] = resolveIdx(k0, normals.size());
           t.n[1] = resolveIdx(k1, normals.size());
           t.n[2] = resolveIdx(k2, normals.size());
       };


       std::string line;
       while (std::getline(f, line)) {
           // Trim trailing whitespace
           line.erase(line.find_last_not_of(" \r\t\n") + 1);
           if (line.empty()) continue;


           const char* s = line.c_str();


           if (line.rfind("usemtl ", 0) == 0) {
               std::string matname = line.substr(7);
               auto result = mtls.emplace(matname, maxGroup + 1);
               if (result.second) {
                   curGroup = ++maxGroup;
               } else {
                   curGroup = result.first->second;
               }
           } else if (line.rfind("vn ", 0) == 0) {
               Vector v;
               sscanf(s, "vn %lf %lf %lf", &v[0], &v[1], &v[2]);
               normals.push_back(v);
           } else if (line.rfind("vt ", 0) == 0) {
               Vector v;
               sscanf(s, "vt %lf %lf", &v[0], &v[1]);
               uvs.push_back(v);
           } else if (line.rfind("v ", 0) == 0) {
               Vector pos, col;
               if (sscanf(s, "v %lf %lf %lf %lf %lf %lf", &pos[0], &pos[1], &pos[2], &col[0], &col[1], &col[2]) == 6) {
                   for (int i = 0; i < 3; i++) col[i] = std::min(1.0, std::max(0.0, col[i]));
                   vertexcolors.push_back(col);
               } else {
                   sscanf(s, "v %lf %lf %lf", &pos[0], &pos[1], &pos[2]);
               }
               vertices.push_back(pos);
           }
           else if (line[0] == 'f') {
               int i[4], j[4], k[4], offset, nn;
               const char* cur = s + 1;
               TriangleIndices t;
               t.group = curGroup;


               // Try each face format: v/vt/vn, v/vt, v//vn, v
               if ((nn = sscanf(cur, "%d/%d/%d %d/%d/%d %d/%d/%d%n", &i[0], &j[0], &k[0], &i[1], &j[1], &k[1], &i[2], &j[2], &k[2], &offset)) == 9) {
                   setFaceVerts(t, i[0], i[1], i[2]);
                   setFaceUVs(t, j[0], j[1], j[2]);
                   setFaceNormals(t, k[0], k[1], k[2]);
               } else if ((nn = sscanf(cur, "%d/%d %d/%d %d/%d%n", &i[0], &j[0], &i[1], &j[1], &i[2], &j[2], &offset)) == 6) {
                   setFaceVerts(t, i[0], i[1], i[2]);
                   setFaceUVs(t, j[0], j[1], j[2]);
               } else if ((nn = sscanf(cur, "%d//%d %d//%d %d//%d%n", &i[0], &k[0], &i[1], &k[1], &i[2], &k[2], &offset)) == 6) {
                   setFaceVerts(t, i[0], i[1], i[2]);
                   setFaceNormals(t, k[0], k[1], k[2]);
               } else if ((nn = sscanf(cur, "%d %d %d%n", &i[0], &i[1], &i[2], &offset)) == 3) {
                   setFaceVerts(t, i[0], i[1], i[2]);
               }
               else continue;


               indices.push_back(t);
               cur += offset;


               // Fan triangulation for polygon faces (4+ vertices)
               while (*cur && *cur != '\n') {
                   TriangleIndices t2;
                   t2.group = curGroup;
                   if ((nn = sscanf(cur, " %d/%d/%d%n", &i[3], &j[3], &k[3], &offset)) == 3) {
                       setFaceVerts(t2, i[0], i[2], i[3]);
                       setFaceUVs(t2, j[0], j[2], j[3]);
                       setFaceNormals(t2, k[0], k[2], k[3]);
                   } else if ((nn = sscanf(cur, " %d/%d%n", &i[3], &j[3], &offset)) == 2) {
                       setFaceVerts(t2, i[0], i[2], i[3]);
                       setFaceUVs(t2, j[0], j[2], j[3]);
                   } else if ((nn = sscanf(cur, " %d//%d%n", &i[3], &k[3], &offset)) == 2) {
                       setFaceVerts(t2, i[0], i[2], i[3]);
                       setFaceNormals(t2, k[0], k[2], k[3]);
                   } else if ((nn = sscanf(cur, " %d%n", &i[3], &offset)) == 1) {
                       setFaceVerts(t2, i[0], i[2], i[3]);
                   } else {
                       cur++;
                       continue;
                   }


                   indices.push_back(t2);
                   cur += offset;
                   i[2] = i[3]; j[2] = j[3]; k[2] = k[3];
               }
           }
       }
   }
  

// empty commit to begin lab 4

 // lab 4 : construction recursive du BVH
 int construireBVH_recursif(int debut, int fin) {
     int idx_noeud = bvh_boites_min.size();

     // calculer la boite englobante de tous les triangles [debut, fin)
     Vector boite_min_locale(1e18, 1e18, 1e18), boite_max_locale(-1e18, -1e18, -1e18);
     for (int k = debut; k < fin; k++) {
         for (int v = 0; v < 3; v++) {
             const Vector& sommet = vertices[indices[k].vtx[v]];
             for (int axe = 0; axe < 3; axe++) {
                 if (sommet[axe] < boite_min_locale[axe]) boite_min_locale[axe] = sommet[axe];
                 if (sommet[axe] > boite_max_locale[axe]) boite_max_locale[axe] = sommet[axe];
             }
         }
     }

     //  now stocker le nouveau noeud dans les vecteurs paralleles
     bvh_boites_min.push_back(boite_min_locale);
     bvh_boites_max.push_back(boite_max_locale);
     bvh_enfants_gauche.push_back(-1);
     bvh_enfants_droite.push_back(-1);
     bvh_debuts_tri.push_back(debut);
     bvh_fins_tri.push_back(fin);

     // critere d'arret : feuille si peu de triangles
     if (fin - debut < 5) return idx_noeud;

     // trouver l'axe le plus long de la diagonale
     Vector diagonale = boite_max_locale - boite_min_locale;
     int axe_plus_long = 0;
     if (diagonale[1] > diagonale[0]) axe_plus_long = 1;
     if (diagonale[2] > diagonale[axe_plus_long]) axe_plus_long = 2;

     double milieu = (boite_min_locale[axe_plus_long] + boite_max_locale[axe_plus_long]) / 2.0;

     // partitionner all triangles autour du milieu 
     int pivot = debut;
     for (int k = debut; k < fin; k++) {
         Vector barycentre = (vertices[indices[k].vtx[0]] + vertices[indices[k].vtx[1]] + vertices[indices[k].vtx[2]]) / 3.0;
         if (barycentre[axe_plus_long] < milieu) {
             TriangleIndices temp = indices[k];
             indices[k] = indices[pivot];
             indices[pivot] = temp;
             pivot++;
         }
     }

     // eviter partition degeneree (tous d'un seul cote)
     if (pivot <= debut || pivot >= fin) pivot = (debut + fin) / 2;

     int idx_gauche = construireBVH_recursif(debut, pivot);
     int idx_droite = construireBVH_recursif(pivot, fin);

     // acceder par index (pas reference) car push_back peut reallouer
     bvh_enfants_gauche[idx_noeud] = idx_gauche;
     bvh_enfants_droite[idx_noeud] = idx_droite;

     return idx_noeud;
 }

 void construireBVH() {
     bvh_boites_min.clear();  bvh_boites_max.clear();
     bvh_enfants_gauche.clear();  bvh_enfants_droite.clear();
     bvh_debuts_tri.clear();  bvh_fins_tri.clear();
     construireBVH_recursif(0, indices.size());
 }

 bool intersect_bvh(const Ray& ray, Vector& P, double& t, Vector& N, int idx_noeud) const {

     // lab 3 : check against the mesh bounding box first (now pernode)
     double tmin_boite = -1e18, tmax_boite = 1e18;
     for (int axe = 0; axe < 3; axe++) {
         double t0 = (bvh_boites_min[idx_noeud][axe] - ray.O[axe]) / ray.u[axe];
         double t1 = (bvh_boites_max[idx_noeud][axe] - ray.O[axe]) / ray.u[axe];
         if (t0 > t1) { double tmp = t0; t0 = t1; t1 = tmp; }
         if (t0 > tmin_boite) tmin_boite = t0;
         if (t1 < tmax_boite) tmax_boite = t1;
     }
     if (tmax_boite < 0 || tmin_boite > tmax_boite) return false;

     // feuille : lab 3 : for each triangle, compute le ray triangle intersection with MollerTrumbore algorithm
     if (bvh_enfants_gauche[idx_noeud] == -1) {
         bool trouve = false;
         for (int k = bvh_debuts_tri[idx_noeud]; k < bvh_fins_tri[idx_noeud]; k++) {
             const Vector& A = vertices[indices[k].vtx[0]];
             const Vector& B = vertices[indices[k].vtx[1]];
             const Vector& C = vertices[indices[k].vtx[2]];

             Vector e1 = B - A;
             Vector e2 = C - A;
             Vector Nvec = cross(e1, e2);

             double denom = dot(ray.u, Nvec);
             if (fabs(denom) < 1e-10) continue;

             Vector AO = A - ray.O;
             Vector AOxu = cross(AO, ray.u);

             double beta  =  dot(e2, AOxu) / denom;
             double gamma = -dot(e1, AOxu) / denom;
             double alpha = 1.0 - beta - gamma;
             double localT = dot(AO, Nvec) / denom;

             if (alpha < 0 || beta < 0 || gamma < 0 || localT < 1e-6) continue;
             if (localT < t) {
                 t = localT;
                 P = ray.O + t * ray.u;
                 N = Nvec;
                 N.normalize();
                 if (dot(ray.u, N) > 0) N = N * (-1.0);
                 trouve = true;
             }
         }
         return trouve;
     }

     // lab 4 : recursively apply the bounding-box test from a BVH datastructure
     bool hit_gauche = intersect_bvh(ray, P, t, N, bvh_enfants_gauche[idx_noeud]);
     bool hit_droite = intersect_bvh(ray, P, t, N, bvh_enfants_droite[idx_noeud]);
     return hit_gauche || hit_droite;
 }


// lab 3 intersect function
 bool intersect(const Ray& ray, Vector& P, double& t, Vector& N) const {
     t = 1e18;
     if (bvh_boites_min.empty()) return false;
     return intersect_bvh(ray, P, t, N, 0);
 }


 std::vector<TriangleIndices> indices;
 std::vector<Vector> vertices;
 std::vector<Vector> normals;
 std::vector<Vector> uvs;
 std::vector<Vector> vertexcolors;
 // lpur lab 4 le BVH stocke dans des vecteurs paralleles (un index= un noeud)
 std::vector<Vector> bvh_boites_min, bvh_boites_max;
 std::vector<int> bvh_enfants_gauche, bvh_enfants_droite;
 std::vector<int> bvh_debuts_tri, bvh_fins_tri;
};

















class Scene {
public:
 Scene() {};
 void addObject(const Object* obj) {
     objects.push_back(obj);
 }








 // returns true iif there is an intersection between the ray and any object in the scene
 // if there is an intersection, also computes the point of the *nearest* intersection P,
 // t>=0 the distance between the ray origin and P (i.e., the parameter along the ray)
 // and the unit normal N.
 // Also returns the index of the object within the std::vector objects in object_id
 bool intersect(const Ray& ray, Vector& P, double& t, Vector& N, int &object_id) const  {








     // TODO (lab 1): iterate through the objects and check the intersections with all of them,
     // and keep the closest intersection, i.e., the one if smallest positive value of t
     t = 1e18;
     bool found = false;
     for (int i = 0; i < (int)objects.size(); i++) {
         Vector localP, localN;
         double localT;
         if (objects[i]->intersect(ray, localP, localT, localN) && localT < t) {
             t = localT;
             P = localP;
             N = localN;
             object_id = i;
             found = true;
         }
     }
     return found;
 }
















 // return the radiance (color) along ray
 Vector getColor(const Ray& ray, int recursion_depth) {








     if (recursion_depth >= max_light_bounce) return Vector(0, 0, 0);








     // TODO (lab 1) : if intersect with ray, use the returned information to compute the color ; otherwise black
     // in lab 1, the color only includes direct lighting with shadows
     Vector P, N;
     double t;
     int object_id;
     if (intersect(ray, P, t, N, object_id)) {








         if (objects[object_id]->mirror) {








             // return getColor in the reflected direction, with recursion_depth+1 (recursively)
         } // else








         if (objects[object_id]->transparent) { // optional








             // return getColor in the refraction direction, with recursion_depth+1 (recursively)
         } // else








         // test if there is a shadow by sending a new ray
         // if there is no shadow, compute the formula with dot products etc.
         Vector light_to_vector = light_position - P;
         Vector light_direction = light_to_vector / light_to_vector.norm();
         Ray shadow_ray(P + 1e-6 * N, light_direction);
         Vector shadow_P, shadow_N;
         double shadow_t;
         int shadow_id;


         Vector direct(0, 0, 0);
         if (intersect(shadow_ray, shadow_P, shadow_t, shadow_N, shadow_id)) {
             if (shadow_t < light_to_vector.norm() - 1e-6) {
                 direct = Vector(0, 0, 0);
             } else {
                 double attenuation = light_intensity / (4 * M_PI * light_to_vector.norm2());
                 Vector material_color = objects[object_id]->albedo / M_PI;
                 double solid_angle = std::max(0., dot(N, light_direction));
                 direct = attenuation * material_color * solid_angle;
             }
         } else {
             double attenuation = light_intensity / (4 * M_PI * light_to_vector.norm2());
             Vector material_color = objects[object_id]->albedo / M_PI;
             double solid_angle = std::max(0., dot(N, light_direction));
             direct = attenuation * material_color * solid_angle;
         }


         // TODO (lab 2) : add indirect lighting component with a recursive call
         int tid = omp_get_thread_num();
         double r1 = uniform(engine[tid]);
         double r2 = uniform(engine[tid]);


         // build local frame around N
         Vector u_dir;
         if (fabs(N[0]) <= fabs(N[1]) && fabs(N[0]) <= fabs(N[2]))
             u_dir = cross(N, Vector(1, 0, 0));
         else if (fabs(N[1]) <= fabs(N[2]))
             u_dir = cross(N, Vector(0, 1, 0));
         else
             u_dir = cross(N, Vector(0, 0, 1));
         u_dir.normalize();
         Vector v_dir = cross(N, u_dir);


         // cosine-weighted random direction in hemisphere
         double cos_t = sqrt(1 - r2);
         double sin_t = sqrt(r2);
         double phi = 2 * M_PI * r1;
         Vector indirect_dir = sin_t * cos(phi) * u_dir + sin_t * sin(phi) * v_dir + cos_t * N;
         indirect_dir.normalize();


         Ray indirect_ray(P + 1e-6 * N, indirect_dir);
         Vector indirect = objects[object_id]->albedo * getColor(indirect_ray, recursion_depth + 1);


         return direct + indirect;
     }








  








     return Vector(0, 0, 0);
 }








 std::vector<const Object*> objects;








 Vector camera_center, light_position;
 double fov, gamma, light_intensity;
 int max_light_bounce;
};
















int main() {
 int W = 512;
 int H = 512;








 for (int i = 0; i<32; i++) {
     engine[i].seed(i);
 }








 Sphere center_sphere(Vector(0, 0, 0), 10., Vector(0.8, 0.8, 0.8));
 Sphere wall_left(Vector(-1000, 0, 0), 940, Vector(0.5, 0.8, 0.1));
 Sphere wall_right(Vector(1000, 0, 0), 940, Vector(0.9, 0.2, 0.3));
 Sphere wall_front(Vector(0, 0, -1000), 940, Vector(0.1, 0.6, 0.7));
 Sphere wall_behind(Vector(0, 0, 1000), 940, Vector(0.8, 0.2, 0.9));
 Sphere ceiling(Vector(0, 1000, 0), 940, Vector(0.3, 0.5, 0.3));
 Sphere floor(Vector(0, -1000, 0), 990, Vector(0.6, 0.5, 0.7));








 Scene scene;
 scene.camera_center = Vector(0, 0, 55);
 scene.light_position = Vector(-10,20,40);
 scene.light_intensity = 3E7;
 scene.fov = 60 * M_PI / 180.;
 scene.gamma = 2.2;    // TODO (lab 1) : play with gamma ; typically, gamma = 2.2
 scene.max_light_bounce = 5;








 TriangleMesh cat(Vector(0.8, 0.6, 0.4));
 cat.readOBJ("cat.obj");
 for (auto& v : cat.vertices) {
     double y = v[1], z = v[2];
     v[1] = z; v[2] = -y;
 }
 for (auto& v : cat.vertices) {
     double x = v[0], z = v[2];
     v[0] = z; v[2] = -x;
 }
 cat.scale_translate(0.6, Vector(0, -10, 0));
 cat.construireBVH(); // lab 4 : construire le BVH apres toutes les transformations
 scene.addObject(&cat);









 scene.addObject(&wall_left);
 scene.addObject(&wall_right);
 scene.addObject(&wall_front);
 scene.addObject(&wall_behind);
 scene.addObject(&ceiling);
 scene.addObject(&floor);








 std::vector<unsigned char> image(W * H * 3, 0);








#pragma omp parallel for schedule(dynamic, 1)
 for (int i = 0; i < H; i++) {
     for (int j = 0; j < W; j++) {
         Vector color;








         // TODO (lab 1) : correct ray_direction so that it goes through each pixel (j, i)        
         Vector ray_direction(j - W/2+0.5, H/2 - i - 0.5, -(W/(2*tan(scene.fov/2))));
         ray_direction.normalize();
         Ray ray(scene.camera_center, ray_direction);








         // TODO (lab 2) : add Monte Carlo / averaging of random ray contributions here
         int tid = omp_get_thread_num();
         int num_mc_trials = 64;
         for (int s = 0; s < num_mc_trials; s++) {


             // TODO (lab 2) : add antialiasing by altering the ray_direction here
             double x_deriv = uniform(engine[tid]) - 0.5;
             double y_deriv = uniform(engine[tid]) - 0.5;
             Vector jittered_direction(j - W/2 + 0.5 + x_deriv, H/2 - i - 0.5 + y_deriv, -(W / (2 * tan(scene.fov / 2))));
             jittered_direction.normalize();
             Ray jittered_ray(scene.camera_center, jittered_direction);


             // TODO (lab 2) : add depth of field effect by altering the ray origin (and direction) here


           color = color + scene.getColor(jittered_ray, 0);
         }
         color = color / num_mc_trials;








         image[(i * W + j) * 3 + 0] = std::min(255., std::max(0., 255. * std::pow(color[0] / 255., 1. / scene.gamma)));
         image[(i * W + j) * 3 + 1] = std::min(255., std::max(0., 255. * std::pow(color[1] / 255., 1. / scene.gamma)));
         image[(i * W + j) * 3 + 2] = std::min(255., std::max(0., 255. * std::pow(color[2] / 255., 1. / scene.gamma)));
     }
 }
 stbi_write_png("image.png", W, H, 3, &image[0], 0);








 return 0;
}





















