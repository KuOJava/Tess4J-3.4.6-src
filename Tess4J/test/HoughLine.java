import java.awt.Graphics2D;
import java.awt.RenderingHints;
import java.awt.image.BufferedImage;
import java.io.FileNotFoundException;
import java.io.FileOutputStream;
import java.io.IOException;
import java.util.ArrayList;

import javax.imageio.ImageIO;


public class HoughLine {

    class h{
        int ro;
        int angle;

        public h(int r,int a){
            this.ro = r;
            this.angle = a;
        }
    }


//    public static void main(String[] args) throws IOException {
//        
//    }

    public BufferedImage hough(BufferedImage im) throws IOException{
        //Image im = this.clone();
        //im.sobel();
        //im.IterBinary();

        int h = im.getHeight();
        int w = im.getWidth();
        int ImageData[] = im.getRGB(0, 0, w, h, null, 0, w);
        int data[] = new int[w*h];

        for(int i=0; i<ImageData.length; i++) {
             data[i] = ImageData[i] & 0xff;
        }

        int ro = (int)Math.sqrt(h*h+w*w);
        int theta = 180;
        int[][] hist = new int[ro][theta];

        for(int k=0;k<theta;k++){
            for(int i=0;i<h;i++){
                for(int j=0;j<w;j++){
                    if(data[j+i*w] != 0){
                        int rho=(int)(j*Math.cos(k*Math.PI/(theta*2))+i*Math.sin(k*Math.PI/(theta*2)));
                        hist[rho][k]++;
                    }      
                }
            }
        }

        ArrayList<h> index = maxIndex(hist,85); //找到大于最大值*0.7的二维直方图的点
        
        double sum=0;


        for(int k =0;k<index.size();k++){

            double resTheta = index.get(k).angle*Math.PI/(theta*2);
            sum+=resTheta;
//            System.out.println(resTheta);

            for(int i=0;i<h;i++){
                for(int j=0;j<w;j++){
                    int rho = (int)(j*Math.cos(resTheta)+i*Math.sin(resTheta));
                    if(data[j+i*w] != 0 && rho == index.get(k).ro){
                        data[j+i*w] = setRed();   //在直线上的点设为红色
                    }else{
                        data[j+i*w] = (255<<24)|(data[j+i*w]<<16)|(data[j+i*w]<<8)|(data[j+i*w]);
                    }   
                }   
            }
        }
        
        sum/=index.size();
        if((Math.PI/2-sum)>0.02)
        	im=rotateImage(im,sum);       
//        im=rotate(im, sum, 0, 0);
        return im;
    }

     public void writeImage(String desImageName, int[] imageData, int width, int height)
     {
      FileOutputStream fos = null;

      try {
       fos = new FileOutputStream(desImageName);
      } catch (FileNotFoundException e) {
       e.printStackTrace();
      }

      BufferedImage img = new BufferedImage(width, height,BufferedImage.TYPE_INT_RGB);
      img.setRGB(0, 0, width, height, imageData, 0, width);
      img.flush();
      try {
       ImageIO.write(img, "jpg", fos);
      } catch (IOException e) {
       e.printStackTrace();
      }
     }

    private ArrayList<HoughLine.h> maxIndex(int[][] hist, int i) {
        ArrayList<h> in = new ArrayList<h>();
        int max = 0;

        for(int i1=0;i1<hist.length;i1++){
            for(int j1=0;j1<hist[i1].length;j1++){
                if(max < hist[i1][j1]){
                    max = hist[i1][j1];
                }
            }
        } 
//        System.out.println(max);

        for(int i1=0;i1<hist.length;i1++){
            for(int j1=0;j1<hist[i1].length;j1++){
                    if(hist[i1][j1] > max*(i/100.0))
                        in.add(new h(i1,j1));
            }
        }

        return in;  
    }

    private int setRed() {

        return 0xffff0000;
    }
    
    public static BufferedImage rotateImage(BufferedImage bufferedimage,  double degree) {  
        int w = bufferedimage.getWidth();  
        int h = bufferedimage.getHeight();  
        int ro = (int)Math.sqrt(h*h+w*w);	
        int type = bufferedimage.getColorModel().getTransparency();  
        BufferedImage img;  
        Graphics2D graphics2d;  
        (graphics2d = (img = new BufferedImage(w, h, type))  
                .createGraphics()).setRenderingHint(  
                RenderingHints.KEY_INTERPOLATION,  
                RenderingHints.VALUE_INTERPOLATION_BILINEAR);  
        graphics2d.rotate(-(Math.PI/2-degree+0.005), w*0.5 , h*0.4 );  
        graphics2d.drawImage(bufferedimage, 0, 0, null);  
        graphics2d.dispose();  
        return img;  
    }
    
    //public static BufferedImage rotate(BufferedImage image, double angle, int cx, int cy) {}

}