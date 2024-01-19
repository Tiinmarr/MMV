#include <iostream>
#include <random>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <filesystem>
#include "qte.h"
#include "implicits.h"
#include "ui_interface.h"
#include "mathematics.h"
#include "mesh.h"
#include "box2.h"
#include "ScalarField.h"
#include "HeightField.h"
#include "perlin.h"


MainWindow::MainWindow() : QMainWindow(), uiw(new Ui::Assets)
{
	// Chargement de l'interface
    uiw->setupUi(this);

	// Chargement du GLWidget
	meshWidget = new MeshWidget;
	QGridLayout* GLlayout = new QGridLayout;
	GLlayout->addWidget(meshWidget, 0, 0);
	GLlayout->setContentsMargins(0, 0, 0, 0);
  uiw->widget_GL->setLayout(GLlayout);

	// Creation des connect
	CreateActions();

	meshWidget->SetCamera(Camera(Vector(10, 0, 0), Vector(0.0, 0.0, 0.0)));
}

MainWindow::~MainWindow()
{
	delete meshWidget;
}

void MainWindow::CreateActions()
{
	// Buttons
	connect(uiw->boxMesh, SIGNAL(clicked()), this, SLOT(BonhommeDeNeige()));
  connect(uiw->resetcameraButton, SIGNAL(clicked()), this, SLOT(ResetCamera()));
  connect(uiw->wireframe, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
  connect(uiw->radioShadingButton_1, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));
  connect(uiw->radioShadingButton_2, SIGNAL(clicked()), this, SLOT(UpdateMaterial()));

	// Widget edition
	connect(meshWidget, SIGNAL(_signalEditSceneLeft(const Ray&)), this, SLOT(editingSceneLeft(const Ray&)));
	connect(meshWidget, SIGNAL(_signalEditSceneRight(const Ray&)), this, SLOT(editingSceneRight(const Ray&)));
}

void MainWindow::editingSceneLeft(const Ray&)
{
}

void MainWindow::editingSceneRight(const Ray&)
{
}

void MainWindow::BonhommeDeNeige()
{

FBM fbm(1.1, 100, 0.65);
FBM fbm2(1.1, 80, 0.62);
FBM fbm3(1.1, 120, 0.7);
RidgeNoise RN(4);

std::random_device rd;
std::mt19937 gen(rd());

Vec2 min(0, 0);
Vec2 max(10, 10);
int rows = 257;
int cols = 257;

HeightField heightField(min, max, rows, cols);
// HeightField noise = heightField;

double centerX1 = 40;
double centerY1 = 40;
double centerX2 = 200;
double centerY2 = 128;
double centerX3 = 128;
double centerY3 = 200;
double centerX4 = 150;
double centerY4 = 60;
double mountainRadius1 = 8.0;
double mountainRadius2 = 20.0;
double mountainRadius3 = 40.0;
double peakRadius = 20.0;


heightField.values[heightField.Index(0,0)] = 5.0;
heightField.values[heightField.Index(0,rows-1)] = 5.0;
heightField.values[heightField.Index(rows-1,0)] = 5.0;
heightField.values[heightField.Index(rows-1,rows-1)] = 5.0;
double h = 100.0;

// Itération sur les niveaux, k = taille d'un carré
for(int k = rows-1; k >= 2; k /= 2, h /= 2.0)
{
int l = k/2; // Demi coté
// Génération pour les carrés
  for(int x=0; x<rows-1;x+=k)
  {
      for(int y=0; y<rows-1; y+=k)
      {
      double a = (heightField.values[heightField.Index(x,y)] +heightField.values[heightField.Index(x+k,y)] +heightField.values[heightField.Index(x,y+k)]+ heightField.values[heightField.Index(x+k,y+k)])/4.0;
      std::uniform_real_distribution<double> distribution(-0.5*h, h);
      double randomNumber = distribution(gen);
      heightField.values[heightField.Index(x+l,y+l)] = a + randomNumber * 0.2;
      }
  }
  for (int x = 0; x < rows - 1; x += l) {
    for (int y = (x+l) % k; y < rows - 1; y += k) {
      double sum =  heightField.values[heightField.Index((x-l+rows)%rows,y)] + heightField.values[heightField.Index((x+l)%rows,y)] +heightField.values[heightField.Index(x,(y+l)%rows)] +  heightField.values[heightField.Index(x,(y-l+rows)%rows)];
      double a = sum / 4.0;
      std::uniform_real_distribution<double> distribution(-0.5*h, h);
      double randomNumber = distribution(gen);
      heightField.values[heightField.Index(x,y)] = a + randomNumber * 0.2;
      if(x == 0) heightField.values[heightField.Index(rows-1,y)] = a;
      if(y == 0) heightField.values[heightField.Index(x,rows-1)] = a;
    }
  }
}

for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
        double distanceToCenter1 = std::sqrt(std::pow(i - centerX1, 2) + std::pow(j - centerY1, 2));
        double distanceToCenter2 = std::sqrt(std::pow(i - centerX2, 2) + std::pow(j - centerY2, 2));
        double distanceToCenter3 = std::sqrt(std::pow(i - centerX3, 2) + std::pow(j - centerY3, 2));
        double distanceToCenter4 = std::sqrt(std::pow(i - centerX4, 2) + std::pow(j - centerY4, 2));
        if (distanceToCenter1 <= mountainRadius1) {
            // À l'intérieur du rayon de la première bosse, générez une hauteur significative
            double heightCenter =  7.0;
            double height = heightCenter + fbm.generate(i, j) - fbm.generate(centerX1, centerY1);
            double slopeFactor = 1.0;
            if (distanceToCenter1 >= 5.0) {
              // Add slope factor
              slopeFactor = 1.0 - ((distanceToCenter1 - 5.0) / (8 - 5.0));
            }
            heightField.values[heightField.Index(i, j)] = heightField.values[heightField.Index(i, j)] + height* slopeFactor;
            
        } else if (distanceToCenter2 <= mountainRadius2) {
            // À l'intérieur du rayon de la deuxième bosse, générez une hauteur plus élevée
            double height = 15.0  + fbm2.generate(i, j) - fbm.generate(centerX2, centerY2);
            double heightFactor = 1.0 - std::exp(-height/ 40.0);
            double slopeFactor = 1.0;
            if (distanceToCenter2 >= 15.0) {
              // Add slope factor
              slopeFactor = 1.0 - ((distanceToCenter2 - 15.0) / (20 - 15.0));
            }
            heightField.values[heightField.Index(i, j)] =heightField.values[heightField.Index(i, j)] + height * heightFactor * slopeFactor;

        } else if (distanceToCenter3 <= peakRadius) {
            // À l'intérieur du rayon de la montagne avec un pic pointu, générez une hauteur maximale
            double height = 50.0 + fbm3.generate(i, j) * 2.0 - fbm.generate(centerX3, centerY3);
            double slopeFactor = 4.0;
            if (distanceToCenter3 >= 5.0) {
              // Add slope factor
              slopeFactor = 1.0 - ((distanceToCenter3 - 5.0) / (20.0 - 5.0));
            }
            heightField.values[heightField.Index(i, j)] = heightField.values[heightField.Index(i, j)] + height * slopeFactor*0.2;
        }
        else if (distanceToCenter4 <= mountainRadius3) {
          // disc primitive : 
          double g = std::pow(1 - (distanceToCenter4 * distanceToCenter4)/ (mountainRadius3 * mountainRadius3),3);
          float e0 = 1 * RN.ridgenoise(i, j);
          float e1 = 0.5 * RN.ridgenoise(2 * i, 2 * j) * e0;
          float e2 = 0.25 * RN.ridgenoise(4 * i, 4 * j) * (e0+e1);
          float e =(e0 + e1 + e2) /(1.75);

          double height;
          if (e<0 ){
            height = 0;
          }
          else {
            height = std::pow(e,3); // elevation
          }
          heightField.values[heightField.Index(i, j)] = g * ( height + 15 ) + heightField.values[heightField.Index(i, j)];
        }
    }
}

// Génération d'un canyon : 
for (int x = 0; x <= 91; x++) {
        int y = 128 - 0.01 * std::pow(x, 2.1);
        if (y < 0) y = 0;
        // Utilisez les coordonnées (x, y) pour intégrer la courbe guide dans votre carte
        heightField.values[heightField.Index(x, y)] -= 25; 
    }
    
    heightField.Smooth();
    // Exporter la carte de hauteurs au format OBJ
    heightField.Export("C:/Users/marti/Desktop/3A_ECL/Master/New MV/MV/Map/HeightField.png",Vector(32,32,40));
    // Carte des pentes :
    heightField.Export_Slope("C:/Users/marti/Desktop/3A_ECL/Master/New MV/MV/Map/Slope.png");
    // Cartes des Laplacien :
    heightField.Export_Laplacian("C:/Users/marti/Desktop/3A_ECL/Master/New MV/MV/Map/Laplacian.png");
    // Cartes de l'accésibilité :
    heightField.Export_Accesibility("C:/Users/marti/Desktop/3A_ECL/Master/New MV/MV/Map/Accessibility.png");

  std::vector<Vector> p;
  std::vector<Vector> nn;
  std::vector<int> t;
 
  p.reserve(rows * rows);
  nn.reserve(rows * rows);
  t.reserve(rows * rows * 6);
  for (int j = 0; j < rows; j++) {
    for (int i = 0; i < rows; i++) {
      p.push_back(heightField.Vertex(i,j));
      nn.push_back(heightField.Normal(i,j));
    }
  }
 
  // Array2 a(Box2::Null, nu, nv);
  Grid a(min, max, rows, cols);
  for (int i = 0; i < rows - 1; i++)
  {
    for (int j = 0; j < rows - 1; j++)
    {
      t.push_back(a.Index(i, j));
      t.push_back(a.Index(i + 1, j));
      t.push_back(a.Index(i + 1, j + 1));
      t.push_back(a.Index(i, j));
      t.push_back(a.Index(i + 1, j + 1));
      t.push_back(a.Index(i, j + 1));
    }
  }
    Mesh mesh(p, nn, t, t);
	
    QString folderPath = "C:\\Users\\marti\\Desktop\\3A_ECL\\Master\\New MV\\MV\\Map\\Map_Mesh2.obj";
	QString fileName = "Map_Mesh";
	mesh.SaveObj(folderPath,fileName);

  std::vector<Color> col;
  col.resize(mesh.Vertexes());
  for (size_t i = 0; i < col.size(); i++)
    col[i] = Color(0.941, 0.902, 0.549);

  meshColor = MeshColor(mesh, col, mesh.VertexIndexes());
  UpdateGeometry();
}

void MainWindow::UpdateGeometry()
{
	meshWidget->ClearAll();
	meshWidget->AddMesh("BoxMesh1", meshColor);

    uiw->lineEdit->setText(QString::number(meshColor.Vertexes()));
    uiw->lineEdit_2->setText(QString::number(meshColor.Triangles()));

	UpdateMaterial();
}

void MainWindow::UpdateMaterial()
{
    meshWidget->UseWireframeGlobal(uiw->wireframe->isChecked());

    if (uiw->radioShadingButton_1->isChecked())
		meshWidget->SetMaterialGlobal(MeshMaterial::Normal);
	else
		meshWidget->SetMaterialGlobal(MeshMaterial::Color);
}

void MainWindow::ResetCamera()
{
	meshWidget->SetCamera(Camera(Vector(-10.0), Vector(0.0)));
}
