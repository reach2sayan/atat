#include <GL/glut.h>

#include <fstream>
#include <iostream>
#include <strstream>

#include "arraylist.h"
#include "getvalue.h"
#include "linalg.h"
#include "meshcalc.h"
#include "meshutil.h"

#define PASTETWO(A, B) A##B

#define MEMBER_CALLBACK(RETTYPE, MEMBER, ARGTYPE, ARGNAME) \
  static RETTYPE PASTETWO(static_, MEMBER) ARGTYPE {       \
    return pcurrent_instance->MEMBER ARGNAME;              \
  }                                                        \
  RETTYPE MEMBER ARGTYPE

#define VOID_MEMBER_CALLBACK(MEMBER, ARGTYPE, ARGNAME) \
  static void PASTETWO(static_, MEMBER) ARGTYPE {      \
    pcurrent_instance->MEMBER ARGNAME;                 \
  }                                                    \
  void MEMBER ARGTYPE

template <class T>
T *is_in_list(const LinkedList<T> &list, const T &val) {
  LinkedListIterator<T> i(list);
  for (; i; i++) {
    if (*i == val) break;
  }
  return i;
}

class BoundingHypercube {
 public:
  Array<Real> min, max;

 public:
  BoundingHypercube(void) : min(), max() {}
  void update(const Array<Real> &v) {
    if (min.get_size() == 0) {
      min = v;
      max = v;
    } else {
      for (int i = 0; i < v.get_size(); i++) {
	if (v(i) < min(i)) min(i) = v(i);
	if (v(i) > max(i)) max(i) = v(i);
      }
    }
  }
};

class ViewerData {
 public:
  static ViewerData *pcurrent_instance;
  LinkedList<int> keysdown;
  int firstdown;
  int modif;
  char lastkey;

 public:
  int dim;
  Array<Array<Real> > pts;
  LinkedList<Array<int> > polys;
  Array<Array<Real> > colors;
  Array2d<Real> rotation;
  Array<Real> center;
  Array<Real> cutpt;
  Array<Array<Real> > old_trans_pts;
  int otheraxis;
  Real rotanglemag;
  Real movemag;
  int writefile;
  int fileindex;

 public:
  ViewerData(void)
      : dim(0),
	pts(),
	polys(),
	rotation(),
	center(),
	keysdown(),
	old_trans_pts() {
    pcurrent_instance = this;
    modif = 0;
    lastkey = 0;
    firstdown = 0;
    rotanglemag = M_PI / 32;
    movemag = 0.02;
    writefile = 0;
    fileindex = 0;
  }
  VOID_MEMBER_CALLBACK(process_normal_keys, (unsigned char key, int x, int y),
		       (key, x, y)) {
    lastkey = key;
  }
  VOID_MEMBER_CALLBACK(press_key, (int key, int x, int y), (key, x, y)) {
    if (!is_in_list(keysdown, key)) {
      keysdown << new int(key);
    }
    modif = glutGetModifiers();
    firstdown = 1;
  }
  VOID_MEMBER_CALLBACK(release_key, (int key, int x, int y), (key, x, y)) {
    LinkedListIterator<int> i(keysdown);
    for (; i; i++) {
      if (*i == key) break;
    }
    if (i) {
      delete keysdown.detach(i);
    }
    modif = glutGetModifiers();
    firstdown = 0;
  }
  VOID_MEMBER_CALLBACK(change_size, (int w, int h), (w, h)) {
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    glViewport(0, 0, w, h);
    gluPerspective(45.0f, (Real)w / (Real)max(1, h), 0.1f, 100.0f);
    glMatrixMode(GL_MODELVIEW);
  }
  VOID_MEMBER_CALLBACK(render_scene, (void), ()) {
    int axis1 = 0;
    int axis2 = 1;
    Real rotangle = 0.;
    Real movedist = 0.;

    int num = lastkey - '1';
    if (num >= 3 && num < dim) {
      otheraxis = num;
    }

    if (lastkey == 'f') {
      rotanglemag *= 2.;
      movemag *= 2.;
    }
    if (lastkey == 's') {
      rotanglemag /= 2.;
      movemag /= 2.;
    }
    if (lastkey == 'r') {
      identity_matrix(&(rotation));
      cutpt = center;
    }
    if (lastkey == 'w') {
      writefile = 1 - writefile;
    }

    Real currotanglemag = rotanglemag;

    if (modif & GLUT_ACTIVE_SHIFT) {
      currotanglemag = M_PI / 2.;
    }

    if (modif & GLUT_ACTIVE_ALT) {
      if (is_in_list(keysdown, GLUT_KEY_UP)) {
	movedist = movemag;
      }
      if (is_in_list(keysdown, GLUT_KEY_DOWN)) {
	movedist = -movemag;
      }
    } else {
      if (is_in_list(keysdown, GLUT_KEY_LEFT) ||
	  is_in_list(keysdown, GLUT_KEY_DOWN) ||
	  is_in_list(keysdown, GLUT_KEY_PAGE_UP)) {
	rotangle = -currotanglemag;
      }
      if (is_in_list(keysdown, GLUT_KEY_RIGHT) ||
	  is_in_list(keysdown, GLUT_KEY_UP) ||
	  is_in_list(keysdown, GLUT_KEY_PAGE_DOWN)) {
	rotangle = currotanglemag;
      }
      if (is_in_list(keysdown, GLUT_KEY_LEFT) ||
	  is_in_list(keysdown, GLUT_KEY_RIGHT)) {
	if (modif & GLUT_ACTIVE_CTRL) {
	  axis1 = 0;
	  axis2 = otheraxis;
	} else {
	  axis1 = 0;
	  axis2 = 2;
	}
      }
      if (is_in_list(keysdown, GLUT_KEY_UP) ||
	  is_in_list(keysdown, GLUT_KEY_DOWN)) {
	if (modif & GLUT_ACTIVE_CTRL) {
	  axis1 = 1;
	  axis2 = otheraxis;
	} else {
	  axis1 = 1;
	  axis2 = 2;
	}
      }
      if (is_in_list(keysdown, GLUT_KEY_PAGE_UP) ||
	  is_in_list(keysdown, GLUT_KEY_PAGE_DOWN)) {
	if (modif & GLUT_ACTIVE_CTRL) {
	  axis1 = 2;
	  axis2 = otheraxis;
	} else {
	  axis1 = 0;
	  axis2 = 1;
	}
      }
    }

    if (fabs(rotangle) > zero_tolerance) {
      //      cout << axis1 << " " << axis2 << " " << rotangle << endl;
    }

    Array2d<Real> drotation(dim, dim);
    identity_matrix(&drotation);
    if (!(modif & GLUT_ACTIVE_SHIFT) || firstdown == 1) {
      drotation(axis1, axis1) = cos(rotangle);
      drotation(axis2, axis2) = cos(rotangle);
      drotation(axis1, axis2) = -sin(rotangle);
      drotation(axis2, axis1) = sin(rotangle);

      cutpt(otheraxis) += movedist;
    }

    Array2d<Real> tmprot(rotation);
    product(&rotation, drotation, tmprot);

    Array<Real> shift;
    product(&shift, rotation, cutpt);
    product(&shift, shift, -1.);
    sum(&shift, shift, cutpt);

    Array<Array<Real> > trans_pts(pts.get_size());
    for (int i = 0; i < trans_pts.get_size(); i++) {
      product(&(trans_pts(i)), rotation, pts(i));
      sum(&trans_pts(i), trans_pts(i), shift);
    }

    Array<Array<Real> > cutpts;
    LinkedList<Array<int> > cutpolys;
    Array<Array<Real> > cutcolors;
    /*
    Array<Array<Real> > cutaxes(3);
    for (int i=0; i<cutaxes.get_size(); i++) {
      cutaxes(i).resize(dim);
      zero_array(&cutaxes(i));
      cutaxes(i)(i)=1.;
    }
    cerr << "begin cut" << endl;
    cross_section_simplexes(&cutpts,&cutpolys,&cutcolors,
    cutaxes,cutpt,trans_pts,polys,colors);
    */
    Array<int> cutaxes(3);
    for (int i = 0; i < cutaxes.get_size(); i++) {
      cutaxes(i) = i;
    }
    // cerr << "begin cut" << endl;
    simple_cross_section_simplexes(&cutpts, &cutpolys, &cutcolors, cutaxes,
				   cutpt, trans_pts, polys, colors);

    // cerr << "pts=" << cutpts.get_size() << " poly=" << cutpolys.get_size() <<
    // endl;
    if (lastkey == 'd') {
      /*
      cout << "diff" << endl;
      for (int i=0; i<old_trans_pts.get_size(); i++) {
	Array<Real> tmp;
	diff(&tmp,old_trans_pts(i),trans_pts(i));
	cout << tmp << endl;
      }
      cout << "pre-CS\n";
      cout << cutpt << endl;
      cout << cutaxes << endl;
      cout << trans_pts<< endl;
      cout << polys << endl;
      cout << "CS\n";
      old_trans_pts=trans_pts;
      */
      cout << cutpts << endl << cutpolys << endl;
    }

    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glLoadIdentity();
    gluLookAt(0.0f, 0.0f, 2.0f, 0.0f, 0.0f, 0.0f, 0.0f, 1.0f, 0.0f);

    glEnable(GL_LIGHTING);
    glEnable(GL_LIGHT0);
    GLfloat lightpos[] = {0.0f, 0.0f, 2.0f, 0.0f};
    glLightfv(GL_LIGHT0, GL_POSITION, lightpos);

    glEnable(GL_CULL_FACE);
    glEnable(GL_DEPTH_TEST);

    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    glBegin(GL_TRIANGLES);

    LinkedListIterator<Array<int> > ipoly(cutpolys);
    int icol = 0;
    for (; ipoly; ipoly++, icol++) {
      Array<rVector3d> border(2);
      for (int s = 0; s < 2; s++) {
	for (int j = 0; j < 3; j++) {
	  border(s)(j) = cutpts((*ipoly)(s + 1))(j) - cutpts((*ipoly)(0))(j);
	}
      }
      rVector3d perp;
      perp = border(0) ^ border(1);
      perp.normalize();
      glNormal3d(perp(0), perp(1), perp(2));

      GLfloat curcol[4];
      curcol[3] = 1.;
      int imod = icol % cutcolors.get_size();
      for (int i = 0; i < cutcolors(imod).get_size(); i++) {
	curcol[i] = cutcolors(imod)(i) / 255.;
      }
      glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, curcol);

      for (int j = 0; j < 3; j++) {
	glVertex3f(cutpts((*ipoly)(j))(0), cutpts((*ipoly)(j))(1),
		   cutpts((*ipoly)(j))(2));
      }

      glNormal3d(-perp(0), -perp(1), -perp(2));
      for (int j = 2; j >= 0; j--) {
	glVertex3f(cutpts((*ipoly)(j))(0), cutpts((*ipoly)(j))(1),
		   cutpts((*ipoly)(j))(2));
      }
    }
    glEnd();

    if (writefile) {
      cerr << "Writing frame" << endl;
      int width = glutGet(GLUT_WINDOW_WIDTH);
      int height = glutGet(GLUT_WINDOW_HEIGHT);
      int bufsize = 3 * width * height;
      GLubyte *data = new GLubyte[bufsize];
      if (data) {
	glReadPixels(0, 0, width, height, GL_RGB, GL_UNSIGNED_BYTE, data);
      }
      ostrstream filename;
      filename << "dumpndviewer" << setfill('0') << setw(4) << fileindex
	       << ".pam"
	       << "\0";
      {
	ofstream dumpfile(filename.str());
	dumpfile << "P7" << endl;
	dumpfile << "WIDTH " << width << endl;
	dumpfile << "HEIGHT " << height << endl;
	dumpfile << "DEPTH 3" << endl;
	dumpfile << "MAXVAL 255" << endl;
	dumpfile << "TUPLTYPE RGB" << endl;
	dumpfile << "ENDHDR" << endl;
	dumpfile.write((char *)data, bufsize);
      }
      delete data;
      fileindex++;
      cerr << "Done writing frame" << endl;
    }
    glutSwapBuffers();
    lastkey = 0;
    firstdown = 0;
  }

  void init(const char *filename) {
    {
      LinkedList<Array<Real> > ptslist;
      LinkedList<Array<Real> > colorslist;
      ifstream file(filename);
      if (!file) {
	ERRORQUIT("Unable to open input file");
      }
      while (1) {
	Array<Array<Real> > curpts;
	LinkedList<Array<int> > curpolys;
	Array<Array<Real> > curcolors;
	cerr << "reading\n";
	file >> curpts;
	if (curpts.get_size() == 0) break;
	file >> curpolys;
	file >> curcolors;
	cerr << curpts.get_size() << " " << curpolys.get_size() << endl;
	combine_mesh(&ptslist, &polys, curpts, curpolys);
	cerr << "done combining\n";
	LinkedListIterator<Array<int> > ip(curpolys);
	for (int i = 0; ip; ip++, i++) {
	  colorslist << new Array<Real>(curcolors(i % curcolors.get_size()));
	}
	cerr << "done colors\n";
      }
      LinkedList_to_Array(&pts, ptslist);
      LinkedList_to_Array(&colors, colorslist);
    }

    cerr << "reading done\n";
    BoundingHypercube box;
    for (int i = 0; i < pts.get_size(); i++) {
      box.update(pts(i));
    }
    sum(&(center), box.min, box.max);
    product(&(center), center, 0.5);
    //    cout << center << endl;
    dim = center.get_size();
    cutpt = center;
    rotation.resize(iVector2d(dim, dim));
    identity_matrix(&(rotation));

    otheraxis = min(3, dim - 1);
    cerr << "setup done\n";
    // init GLUT and create window
    glutInitDisplayMode(GLUT_DEPTH | GLUT_DOUBLE | GLUT_RGBA);
    glutInitWindowPosition(100, 100);
    glutInitWindowSize(600, 600);
    glutCreateWindow("ndviewer");

    glutIgnoreKeyRepeat(1);
    glutDisplayFunc(&ViewerData::static_render_scene);
    glutIdleFunc(&ViewerData::static_render_scene);
    glutReshapeFunc(&ViewerData::static_change_size);
    glutKeyboardFunc(&ViewerData::static_process_normal_keys);
    glutSpecialFunc(&ViewerData::static_press_key);
    glutSpecialUpFunc(&ViewerData::static_release_key);
    cerr << "glut init done\n";
  }
};

ViewerData *ViewerData::pcurrent_instance = NULL;

int main(int argc, char *argv[]) {
  char *infilename = "mesh.nd";
  AskStruct options[] = {
      {"", "N-Dimensional VIEWER " MAPS_VERSION ", by Axel van de Walle",
       TITLEVAL, NULL},
      {"-if", "Input table file", STRINGVAL, &infilename}};
  if (!get_values(argc, argv, countof(options), options)) {
    display_help(countof(options), options);
    return 1;
  }
  /*
  { //debug;
    Array<Array<Real> > pts;
    ifstream file("pt.in");
    file >> pts;
    Array<Array<Real> > pts3(3);
    pts3(0)=pts(0);
    pts3(1)=pts(1);
    pts3(2)=pts(3);
    Array<Real> center;
    Real radius2;
    calc_circumsphere(&center, &radius2,pts3);
    cout << center << endl;
    cout << radius2 << endl;
    LinkedList<Array<int> > polys;
    create_mesh(&polys, pts,0,MAXFLOAT);
    cout << polys;
    return 1;
  }
  */
  crash_if_singular = 0;
  glutInit(&argc, argv);
  ViewerData viewerdata;
  viewerdata.init(infilename);
  cerr << "Starting main loop\n";
  glutMainLoop();

  /*
  Array<Real> cutpt;
  Array<Array<Real> > cutaxes;
  Array<Array<Real> > pts;
  LinkedList<Array<int> > polys;
  cin >> cutpt;
  cin >> cutaxes;
  cin >> pts;
  cin >> polys;
  Array<Array<Real> > cutpts;
  LinkedList<Array<int> > cutpolys;
  cross_section_simplexes(&cutpts,&cutpolys, cutaxes,cutpt,pts,polys);
  cout << cutpts << endl;
  cout << cutpolys << endl;
  return 1;
  */

  return 0;
}
