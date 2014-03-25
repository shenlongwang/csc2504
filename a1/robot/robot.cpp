/***********************************************************
             CSC418/2504, Fall 2009
  
                 penguin.cpp
                 
       Simple demo program using OpenGL and the glut/glui 
       libraries

  
    Instructions:
        Please read the assignment page to determine 
        exactly what needs to be implemented.  Then read 
        over this file and become acquainted with its 
        design.

        Add source code where it appears appropriate. In
        particular, see lines marked 'TODO'.

        You should not need to change the overall structure
        of the program. However it should be clear what
        your changes do, and you should use sufficient comments
        to explain your code.  While the point of the assignment
        is to draw and animate the character, you will
        also be marked based on your design.

***********************************************************/

#ifdef _WIN32
#include <windows.h>
#endif

#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>
#include <GL/glui.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifndef _WIN32
#include <unistd.h>
#else
void usleep(unsigned int nanosec)
{
    Sleep(nanosec / 1000);
}
#endif


// *************** GLOBAL VARIABLES *************************


const float PI = 3.14159;

// --------------- USER INTERFACE VARIABLES -----------------

// Window settings
int windowID;               // Glut window ID (for display)
GLUI *glui;                 // Glui window (for controls)
int Win[2];                 // window (x,y) size


// ---------------- ANIMATION VARIABLES ---------------------

// Animation settings
int animate_mode = 0;       // 0 = no anim, 1 = animate
int animation_frame = 0;      // Specify current frame of animation

// Joint parameters
const float JOINT_MIN = -45.0f;
const float JOINT_MAX =  45.0f;
float joint_rot = 0.0f;

//////////////////////////////////////////////////////
// TODO: Add additional joint parameters here
//////////////////////////////////////////////////////
// Declare parameters: 

// Minimum and Maximum of Arm/Leg Length
const float LENGTH_MIN = 20.0f;
const float LENGTH_MAX = 60.0f;
// Initial rotation angles of left arm/ left leg/ left upperarm
float leftarm_rot = 0.0f;
float leftleg_rot = 0.0f;
float leftupperarm_rot = 0.0f;
// Initial rotation angles of right arm/ right leg/ right upperarm
float rightarm_rot = 0.0f;
float rightleg_rot = 0.0f;
float rightupperarm_rot = 0.0f;
// Initial length  of right upperarm/ right leg/left upperarm/ left leg 
float leftleg_len = 20.0f;
float leftupperarm_len = 20.0f;
float rightleg_len = 20.0f;
float rightupperarm_len = 20.0f;

// ***********  FUNCTION HEADER DECLARATIONS ****************


// Initialization functions
void initGlut(char* winName);
void initGlui();
void initGl();


// Callbacks for handling events in glut
void myReshape(int w, int h);
void animate();
void display(void);

// Callback for handling events in glui
void GLUI_Control(int id);


// Functions to help draw the object
void drawSquare(float size);
void drawCircle(float radius);


// Return the current system clock (in seconds)
double getTime();


// ******************** FUNCTIONS ************************



// main() function
// Initializes the user interface (and any user variables)
// then hands over control to the event handler, which calls 
// display() whenever the GL window needs to be redrawn.
int main(int argc, char** argv)
{

    // Process program arguments
    if(argc != 3) {
        printf("Usage: demo [width] [height]\n");
        printf("Using 300x200 window by default...\n");
        Win[0] = 300;
        Win[1] = 200;
    } else {
        Win[0] = atoi(argv[1]);
        Win[1] = atoi(argv[2]);
    }


    // Initialize glut, glui, and opengl
    glutInit(&argc, argv);
    initGlut(argv[0]);
    initGlui();
    initGl();

    // Invoke the standard GLUT main event loop
    glutMainLoop();

    return 0;         // never reached
}


// Initialize glut and create a window with the specified caption 
void initGlut(char* winName)
{
    // Set video mode: double-buffered, color, depth-buffered
    glutInitDisplayMode (GLUT_DOUBLE | GLUT_RGB | GLUT_DEPTH);

    // Create window
    glutInitWindowPosition (0, 0);
    glutInitWindowSize(Win[0],Win[1]);
    windowID = glutCreateWindow(winName);

    // Setup callback functions to handle events
    glutReshapeFunc(myReshape); // Call myReshape whenever window resized
    glutDisplayFunc(display);   // Call display whenever new frame needed 
}


// Quit button handler.  Called when the "quit" button is pressed.
void quitButton(int)
{
  exit(0);
}

// Animate button handler.  Called when the "animate" checkbox is pressed.
void animateButton(int)
{
  // synchronize variables that GLUT uses
  glui->sync_live();

  animation_frame = 0;
  if(animate_mode == 1) {
    // start animation
    GLUI_Master.set_glutIdleFunc(animate);
  } else {
    // stop animation
    GLUI_Master.set_glutIdleFunc(NULL);
  }
}

// Initialize GLUI and the user interface
void initGlui()
{
    GLUI_Master.set_glutIdleFunc(NULL);

    // Create GLUI window
    glui = GLUI_Master.create_glui("Glui Window", 0, Win[0]+10, 0);

    // Create a control to specify the rotation of the joint
    GLUI_Spinner *joint_spinner
        = glui->add_spinner("Neck", GLUI_SPINNER_FLOAT, &joint_rot);
    joint_spinner->set_speed(0.1);
    joint_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    ///////////////////////////////////////////////////////////
    // TODO: 
    //   Add controls for additional joints here
    ///////////////////////////////////////////////////////////

    // Create a control to specify the rotation of the left arm 
    GLUI_Spinner *leftarmrot_spinner
        = glui->add_spinner("Left Arm Rotation", GLUI_SPINNER_FLOAT, &leftarm_rot);
    leftarmrot_spinner->set_speed(0.1);
    leftarmrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the left leg 
    GLUI_Spinner *leftlegrot_spinner
        = glui->add_spinner("Left Leg Rotation", GLUI_SPINNER_FLOAT, &leftleg_rot);
    leftlegrot_spinner->set_speed(0.1);
    leftlegrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the left upper arm 
    GLUI_Spinner *leftupperrot_spinner
        = glui->add_spinner("Left Upper Arm Rotation", GLUI_SPINNER_FLOAT, &leftupperarm_rot);
    leftupperrot_spinner->set_speed(0.1);
    leftupperrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the right arm 
    GLUI_Spinner *rightarmrot_spinner
        = glui->add_spinner("Right Arm Rotation", GLUI_SPINNER_FLOAT, &rightarm_rot);
    rightarmrot_spinner->set_speed(0.1);
    rightarmrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the right leg 
    GLUI_Spinner *rightlegrot_spinner
        = glui->add_spinner("Right Leg Rotation", GLUI_SPINNER_FLOAT, &rightleg_rot);
    rightlegrot_spinner->set_speed(0.1);
    rightlegrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the right upper arm 
    GLUI_Spinner *rightupperarmrot_spinner
        = glui->add_spinner("Right Upper Arm Rotation", GLUI_SPINNER_FLOAT, &rightupperarm_rot);
    rightupperarmrot_spinner->set_speed(0.1);
    rightupperarmrot_spinner->set_float_limits(JOINT_MIN, JOINT_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the left upper arm 
    GLUI_Spinner *leftupperarmlen_spinner
        = glui->add_spinner("Left Upper Arm Length", GLUI_SPINNER_FLOAT, &leftupperarm_len);
    leftupperarmlen_spinner->set_speed(0.1);
    leftupperarmlen_spinner->set_float_limits(LENGTH_MIN, LENGTH_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the left leg 
    GLUI_Spinner *leftleglen_spinner
        = glui->add_spinner("Left Leg Length", GLUI_SPINNER_FLOAT, &leftleg_len);
    leftleglen_spinner->set_speed(0.1);
    leftleglen_spinner->set_float_limits(LENGTH_MIN, LENGTH_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the right upper arm 
    GLUI_Spinner *rightupperarmlen_spinner
        = glui->add_spinner("Right Upper Arm Length", GLUI_SPINNER_FLOAT, &rightupperarm_len);
    rightupperarmlen_spinner->set_speed(0.1);
    rightupperarmlen_spinner->set_float_limits(LENGTH_MIN, LENGTH_MAX, GLUI_LIMIT_CLAMP);

    // Create a control to specify the rotation of the Right leg
    GLUI_Spinner *rightleglen_spinner
        = glui->add_spinner("Right Leg Length", GLUI_SPINNER_FLOAT, &rightleg_len);
    rightleglen_spinner->set_speed(0.1);
    rightleglen_spinner->set_float_limits(LENGTH_MIN, LENGTH_MAX, GLUI_LIMIT_CLAMP);

    // Add button to specify animation mode 
    glui->add_separator();
    glui->add_checkbox("Animate", &animate_mode, 0, animateButton);

    // Add "Quit" button
    glui->add_separator();
    glui->add_button("Quit", 0, quitButton);

    // Set the main window to be the "active" window
    glui->set_main_gfx_window(windowID);
}


// Performs most of the OpenGL intialization
void initGl(void)
{
    // glClearColor (red, green, blue, alpha)
    // Ignore the meaning of the 'alpha' value for now
    glClearColor(0.7f,0.7f,0.9f,1.0f);
}




// Callback idle function for animating the scene
void animate()
{
    // Update geometry
    const double joint_rot_speed = 0.1;
    double joint_rot_t = (sin(animation_frame*joint_rot_speed) + 1.0) / 2.0;
    joint_rot = joint_rot_t * JOINT_MIN + (1 - joint_rot_t) * JOINT_MAX;
    
    ///////////////////////////////////////////////////////////
    // TODO:
    //   Modify this function animate the character's joints
    //   Note: Nothing should be drawn in this function!  OpenGL drawing
    //   should only happen in the display() callback.
    ///////////////////////////////////////////////////////////
    leftarm_rot = joint_rot;
    leftleg_rot = joint_rot;
    leftupperarm_rot = joint_rot;
    rightarm_rot = joint_rot;
    rightleg_rot = joint_rot;
    rightupperarm_rot = joint_rot;
    leftleg_len = abs(joint_rot)*PI/180*50 + 20.0;
    leftupperarm_len = abs(joint_rot)*PI/180*50 + 20.0;
    rightleg_len = abs(joint_rot)*PI/180*50 + 20.0;
    rightupperarm_len = abs(joint_rot)*PI/180*50 + 20.0;

    // Update user interface
    glui->sync_live();

    // Tell glut window to update itself.  This will cause the display()
    // callback to be called, which renders the object (once you've written
    // the callback).
    glutSetWindow(windowID);
    glutPostRedisplay();

    // increment the frame number.
    animation_frame++;

    // Wait 50 ms between frames (20 frames per second)
    usleep(50000);
}


// Handles the window being resized by updating the viewport
// and projection matrices
void myReshape(int w, int h)
{
    // Setup projection matrix for new window
    glMatrixMode(GL_PROJECTION);
    glLoadIdentity();
    gluOrtho2D(-w/2, w/2, -h/2, h/2);

    // Update OpenGL viewport and internal variables
    glViewport(0,0, w,h);
    Win[0] = w;
    Win[1] = h;
}



// display callback
//
// This gets called by the event handler to draw
// the scene, so this is where you need to build
// your scene -- make your changes and additions here.
// All rendering happens in this function.  For Assignment 1,
// updates to geometry should happen in the "animate" function.
void display(void)
{
    // glClearColor (red, green, blue, alpha)
    // Ignore the meaning of the 'alpha' value for now
    glClearColor(0.7f,0.7f,0.9f,1.0f);

    // OK, now clear the screen with the background colour
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);  

    // Setup the model-view transformation matrix
    glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();

    ///////////////////////////////////////////////////////////
    // TODO:
    //   Modify this function draw the scene
    //   This should include function calls to pieces that
    //   apply the appropriate transformation matrice and
    //   render the individual body parts.
    ///////////////////////////////////////////////////////////

    // Draw our hinged object
    const float BODY_WIDTH = 30.0f;
    const float BODY_LENGTH = 50.0f;
    const float ARM_LENGTH = 50.0f;
    const float ARM_WIDTH = 10.0f;
    // Push the current transformation matrix on the stack
    glPushMatrix();

        // Draw the 'body'
        // Scale square to size of body
   
        glScalef(BODY_WIDTH, BODY_LENGTH, 1.0);
        // Set the colour to green
        glColor3f(0.0, 1.0, 0.0);
        drawSquare(1.0);
 
        // Draw the 'lower body'
        
	glPushMatrix(); 
            // Move the lower body to the hinge
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
	    glTranslatef(0.0, -BODY_LENGTH/2 -ARM_WIDTH/2, 0.0);
	    // Scale square to size of body
            glScalef(ARM_LENGTH, ARM_WIDTH, 1.0);

            // Set the colour to green
            glColor3f(0.0, 1.0, 1.0);

            // Draw the square for the body
            drawSquare(1.0);
       
            // Draw the 'right leg'
    	    glPushMatrix(); 
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the hinge
                glTranslatef(-BODY_WIDTH/2 - ARM_WIDTH/2, 0.0, 0.0);
                // Rotate along the hinge
                glRotatef(rightleg_rot, 0.0, 0.0, 1.0);
        	// Move to center location of arm, under previous rotation
                glTranslatef(0.0, -rightleg_len/2+ARM_WIDTH/2, 0.0);
                // Scale the size of the arm
                glScalef(ARM_WIDTH, rightleg_len, 1.0);
                // Set the colour to red 
                glColor3f(1.0, 0.0, 0.0);
                // Draw the square for the leg 
                drawSquare(1.0);
            glPopMatrix();
            
	    // Draw the 'rightleg circle'
    	    glPushMatrix(); 
                // Draw the 'rightleg'
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the joint hinge
                glTranslatef(-BODY_WIDTH/2 - ARM_WIDTH/2, 0.0, 0.0);
                // Set the colour to black 
                glColor3f(0.0, 0.0, 0.0);
                // Draw the circle 
                drawCircle(2.0); 
            glPopMatrix();
            
	    // Draw the 'left leg'
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
		// Move the arm to the hinge
                glTranslatef(BODY_WIDTH/2 + ARM_WIDTH/2, 0.0, 0.0);
                // Rotate along the hinge
                glRotatef(leftleg_rot, 0.0, 0.0, -1.0);
        	// Move to center location of arm, under previous rotation
                glTranslatef(0.0, -leftleg_len/2 + ARM_WIDTH/2, 0.0);
		// Scale the size of the arm
                glScalef(ARM_WIDTH, leftleg_len, 1.0);
                // Set the colour to red 
                glColor3f(1.0, 0.0, 0.0);
                // Draw the square for the leg 
                drawSquare(1.0);
    	    glPopMatrix();
	  
            // Draw the 'left leg circle'
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
		// Move the arm to the joint hinge
                glTranslatef(BODY_WIDTH/2 + ARM_WIDTH/2, 0.0, 0.0);
                // Set the colour to black 
                glColor3f(0.0, 0.0, 0.0);
                // Draw the circle 
                drawCircle(2.0); 
    	    glPopMatrix();
        glPopMatrix();
      
        // Draw the 'leftarm'
        glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
	    // Move the arm to the joint hinge
            glTranslatef(-BODY_WIDTH/2, ARM_WIDTH/2, 0.0);
            // Rotate along the hinge
            glRotatef(leftarm_rot, 0.0, 0.0, -1.0);
            // Move to center location of arm, under previous rotation
            glTranslatef(-ARM_LENGTH/2+ARM_WIDTH/2, 0.0, 0.0);
            // Scale the size of the arm
            glScalef(ARM_LENGTH, ARM_WIDTH, 1.0);
            // Set the colour to red 
            glColor3f(1.0, 0.0, 0.0);
            // Draw the square for the arm
            drawSquare(1.0);
    
   	    // Draw left upper arm	 
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the joint hinge
                glTranslatef(-ARM_LENGTH/2+ARM_WIDTH/2, 0.0, 0.0);
                // Rotate along the hinge
                glRotatef(leftupperarm_rot, 0.0, 0.0, -1.0);
                // Move to center location of arm, under previous rotation
                glTranslatef(-leftupperarm_len/2+ARM_WIDTH/2, 0.0, 0.0);
                // Scale the size of the arm
                glScalef(leftupperarm_len, ARM_WIDTH, 1.0);
                // set the colour to black 
                glColor3f(1.0, 0.0, 1.0);
                // Draw the square for the arm
                drawSquare(1.0); 
            glPopMatrix();
    
   	    // Draw Upper Arm Circle	 
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the joint hinge
                glTranslatef(-ARM_LENGTH/2+ARM_WIDTH/2, 0.0, 0.0);
                // set the colour to black 
                glColor3f(0.0, 0.0, 0.0);
                // Draw the circle 
                drawCircle(2.0); 
            glPopMatrix();
    	glPopMatrix();
 
        // Draw the 'leftarm circle'
        glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
    	    // move the arm to the joint hinge
            glTranslatef(-BODY_WIDTH/2, ARM_WIDTH/2, 0.0);
            // Set the colour to black 
            glColor3f(0.0, 0.0, 0.0);
            // Draw the circle 
            drawCircle(2.0); 
        glPopMatrix();
        // Draw the 'rightarm'
        glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
	    // Move the arm to the joint hinge
            glTranslatef(BODY_WIDTH/2, ARM_WIDTH/2, 0.0);
            // Rotate along the hinge
            glRotatef(rightarm_rot, 0.0, 0.0, -1.0);
    	    // Move to center location of arm, under previous rotation
            glTranslatef(ARM_LENGTH/2-ARM_WIDTH/2, 0.0, 0.0);
            // Scale the size of the arm
            glScalef(ARM_LENGTH, ARM_WIDTH, 1.0);
            // Draw the square for the arm
            glColor3f(1.0, 0.0, 0.0);
            // Draw the square for the arm
            drawSquare(1.0);
    
   	    // Draw right upper arm	 
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the joint hinge
                glTranslatef(ARM_LENGTH/2-ARM_WIDTH/2, 0.0, 0.0);
                // Rotate along the hinge
                glRotatef(rightupperarm_rot, 0.0, 0.0, -1.0);
                // Move to center location of arm, under previous rotation
                glTranslatef(rightupperarm_len/2-ARM_WIDTH/2, 0.0, 0.0);
                // Scale the size of the arm
                glScalef(rightupperarm_len, ARM_WIDTH, 1.0);
                // Set the colour to yellow 
                glColor3f(1.0, 1.0, 0.0);
                // Draw the square for the arm
                drawSquare(1.0); 
            glPopMatrix();
    
            glPushMatrix();
                // Scale them back
                glScalef(1/ARM_LENGTH, 1/ARM_WIDTH, 1.0);
        	// Move the arm to the joint hinge
                glTranslatef(ARM_LENGTH/2-ARM_WIDTH/2, 0.0, 0.0);
                // Set the colour to black 
                glColor3f(0.0, 0.0, 0.0);
                // Draw the circle 
                drawCircle(2.0); 
            glPopMatrix();
	glPopMatrix();
 
        // Draw the 'rightarm circle'
        glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
    	    // Move the arm to the joint hinge
            glTranslatef(BODY_WIDTH/2, ARM_WIDTH/2, 0.0);
            // Set the colour to black 
            glColor3f(0.0, 0.0, 0.0);
            // Draw the circle 
            drawCircle(2.0); 
        glPopMatrix();
 
	// Draw the 'neck'

        // Move the arm to the joint hinge
        glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
    	    // Move the arm to the joint hinge
            glTranslatef(0.0, BODY_LENGTH/2 - ARM_WIDTH, 0.0);
            // Rotate along the hinge
            glRotatef(joint_rot, 0.0, 0.0, 1.0);
            // Move to center location of arm, under previous rotation
            glTranslatef(0.0, ARM_LENGTH/2 - ARM_WIDTH/2, 0.0);
            // Scale the size of the neck 
            glScalef(ARM_WIDTH, ARM_LENGTH, 1.0);
            // Set the colour to yellow 
            glColor3f(1.0, 1.0, 0.0);
            // Draw the square for the neck 
            drawSquare(1.0);
    
            // Move the arm to the joint hinge
            glPushMatrix();
                // Scale them back
    		glScalef(1/ARM_WIDTH, 1/ARM_LENGTH, 1.0); 
    	        // Move the arm to the joint hinge
		glTranslatef(0.0, BODY_LENGTH/2 - ARM_WIDTH, 0.0);
              	// Scale the size of the head 
                glScalef(ARM_LENGTH, ARM_LENGTH/2, 1.0);
                // Set the colour to yellow 
                glColor3f(1.0, 1.0, 0.0);
                // Draw the square for the head
                drawSquare(1.0);

		// Draw the 'left eye'
                glPushMatrix();
                    // Scale them back
                    glScalef(1/ARM_LENGTH, 2/ARM_LENGTH, 1.0);
    	            // Move the arm to the joint hinge
		    glTranslatef(-15.0, 0.0, 0.0);
                    // Set the colour to yellow 
                    glColor3f(0.0, 0.0, 0.0);
                    // Draw the big circle 
                    drawCircle(10.0);
                    // Draw the small circle 
                    drawCircle(5.0);
                glPopMatrix();	
        
        	// Draw the 'right eye'
                glPushMatrix();
                    // Scale them back
                    glScalef(1/ARM_LENGTH, 2/ARM_LENGTH, 1.0);
    	            // Move the arm to the joint hinge
		    glTranslatef(15.0, 0.0, 0.0);
                    // Set the colour to yellow 
                    glColor3f(0.0, 0.0, 0.0);
                    // Draw the big circle 
                    drawCircle(10.0);
                    // Draw the small circle 
                    drawCircle(5.0);
                glPopMatrix();	
    
            glPopMatrix();	
            
        glPopMatrix();	
        // Draw neck circle
	glPushMatrix();
            // Scale them back
            glScalef(1/BODY_WIDTH, 1/BODY_LENGTH, 1.0);
    	    // Move the arm to the joint hinge
            glTranslatef(0.0, BODY_LENGTH/2 - ARM_WIDTH, 0.0);
            // Set the colour to black 
            glColor3f(0.0, 0.0, 0.0);
            // Draw the circle 
            drawCircle(2.0); 
        glPopMatrix();
    glPopMatrix();
    // Retrieve the previous state of the transformation stack

    // Execute any GL functions that are in the queue just to be safe
    glFlush();

    // Now, show the frame buffer that we just drew into.
    // (this prevents flickering).
    glutSwapBuffers();
}


// Draw a square of the specified size, centered at the current location
void drawSquare(float width)
{
    // Draw the square
    glBegin(GL_POLYGON);
    glVertex2d(-width/2, -width/2);
    glVertex2d(width/2, -width/2);
    glVertex2d(width/2, width/2);
    glVertex2d(-width/2, width/2);
    glEnd();
}

// Draw a circle of the specified radius, centered at the current location
const float DEG2RAD = 3.14159/180;
 
void drawCircle(float radius)
{
   glBegin(GL_LINE_LOOP);
 
   for (int i=0; i<360; i++)
   {
      float degInRad = i*DEG2RAD;
      glVertex2f(cos(degInRad)*radius,sin(degInRad)*radius);
   }
 
   glEnd();
}
