//package LSD;
//
//
//import javax.swing.*;
//import javax.imageio.ImageIO;
//import java.awt.image.BufferedImage;
//import java.io.File;
//import java.io.IOException;
//import java.util.HashSet;
//import java.awt.*;
//import java.awt.image.*;
//import javax.swing.*;
//
//
//
//public class GUI extends JFrame {
//
//
//	GUI() {
//		setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
//		double max_x1=0;
//		double max_x2=0;
//		double max_y1=0;
//		double max_y2=0;
//		double min_x1=0;
//		double min_x2=0;
//		double min_y1=0;
//		double min_y2=0;
//
//		double min_x11=0;
//		double min_x22=0;
//		double min_y11=0;
//		double min_y22=0;
//
//		try {
//			BufferedImage myPicture = ImageIO.read(new File("././b/3.jpg"));
//			Graphics2D g2d = myPicture.createGraphics();
//			int x = myPicture.getWidth();
//			int y = myPicture.getHeight();
//			max_x1=x;
//			System.out.println(x);
//			System.out.println(y);
//		
//			HashSet<Line> lines = new HashSet<Line>();
//
//			
//			double [] arr = myPicture.getData().getPixels(0,0,x,y,new double[x*y*3]);
//
//			double [] arr2 = new double[x*y];
//		
//			System.out.println("arr.length:"+arr.length);
//			int c=0;
//			for(int i = 0; i < arr.length-3; i+=3) {
//				double B = arr[i];
//				double G = arr[i+1];
//				double R = arr[i+2];
//				double level = R * 0.2126 + G * 0.7152 + B * 0.0722;
//				arr2[c++] = level;
//			}
//
//			LSD lsd = new LSD();
//
//			double [] out = lsd.lsd(arr2,x,y);
//
//			for(int i = 0; i < lsd.n_out; i++) {
//				for (int j = 0; j < 7; j++)
//				
//				lines.add(new Line(out[7 * i + 0], out[7 * i + 1],
//						out[7 * i + 2], out[7 * i + 3]));
//
//			}
//			
//			
//			int z=0;
//			System.out.println("¿ªÊ¼");
//			for (Line l : lines) {
//				if((l.x2-l.x1>x/8)&&(l.y1>min_y11)&&(l.y1<y/2)&&(Math.abs(l.y1-l.y2)<10))
//				{
//					System.out.println("true");
//					min_x11=l.x1;
//					min_x22=l.x2;
//					min_y11=l.y1;
//					min_y22=l.y2;
//				}
//				if(l.x2-l.x1>x/8)
//				{	
//				System.out.println(l.y2-l.y1);
//				System.out.println(l.x1);
//				System.out.println(l.x2);
//				System.out.println(l.y1);
//				System.out.println(l.y2);
//				System.out.println(++z);}
//				if((l.y2-l.y1>y/8)&&(l.x1>min_x1)&&(l.x1<x/2)&&(Math.abs(l.x1-l.x2)<10))
//				{
//					min_x1=l.x1;
//					min_x2=l.x2;
//					min_y1=l.y1;
//					min_y2=l.y2;
//				}
//				else if((l.y2-l.y1>y/8)&&(l.x1<max_x1)&&(l.x1>x/2)&&(Math.abs(l.x1-l.x2)<10))
//				{
//					
//					max_x1=l.x1;
//					max_x2=l.x2;
//					max_y1=l.y1;
//					max_y2=l.y2;
//				}
//				
//			}
//			g2d.setPaint(Color.red);
//		//	g2d.drawLine((int)max_x1,(int)max_y1,(int)max_x2,(int)max_y2);
//		
//		
//			//g2d.drawLine((int)min_x1,(int)min_y1,(int)min_x2,(int)min_y2);
//			g2d.drawLine((int)min_x11,(int)min_y11,(int)min_x22,(int)min_y22);
//		
//			System.out.println(max_x1+" "+min_x1+" "+min_y11);
//			System.out.println(x/2);
//			System.out.println(y/8);
//			System.out.println("½áÊø");
//			System.out.println(lines.size());
//			BufferedImage img=myPicture.getSubimage ((int)min_x1,(int)min_y11,(int)(max_x1-min_x1),(int)(y*2/3));
//			JLabel picLabel = new JLabel(new ImageIcon(img));
//			add(picLabel);
//
//
//		} catch (IOException e) {
//
//		}
//
//		setSize(600,600);
//        	setVisible(true);
//	}
//	public static void main(String [] args){
//		new GUI();
//	
//	}
//
//
//}
