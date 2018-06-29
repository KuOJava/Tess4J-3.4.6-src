import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.IOException;

public class Invert {

	public static BufferedImage invert(BufferedImage img) throws IOException{
		int h=img.getHeight();
		int w=img.getWidth();
		int[][] pixels=new int[w][h];
		for(int i=0;i<w;i++){
			for(int j=0;j<h;j++){
				pixels[i][j]=setwhite(img.getRGB(i, j));
			}
		}
		BufferedImage nbi=new BufferedImage(w,h,BufferedImage.TYPE_BYTE_BINARY);
		int sw=-6579300;
		for (int x = 0; x < w; x++) {  
            for (int y = 0; y < h; y++) {  
                if(pixels[x][y]<sw){  
                    int max=new Color(255,255,255).getRGB();  
                    nbi.setRGB(x, y, max);  
                }else{  
                    int min=new Color(0,0,0).getRGB();  
                    nbi.setRGB(x, y, min);  
                }  
            }  
        }    		
        return nbi;
	}
	private static int getAverageColor(int[][] gray, int x, int y, int w, int h) {
		int rs = gray[x][y]  
                + (x == 0 ? 255 : gray[x - 1][y])  
                + (x == 0 || y == 0 ? 255 : gray[x - 1][y - 1])  
                + (x == 0 || y == h - 1 ? 255 : gray[x - 1][y + 1])  
                + (y == 0 ? 255 : gray[x][y - 1])  
                + (y == h - 1 ? 255 : gray[x][y + 1])  
                + (x == w - 1 ? 255 : gray[x + 1][ y])  
                + (x == w - 1 || y == 0 ? 255 : gray[x + 1][y - 1])  
                + (x == w - 1 || y == h - 1 ? 255 : gray[x + 1][y + 1]);  
		return rs / 9;
	}

	private static int getGray(int rgb) {
		Color c=new Color(rgb);
		int r=c.getRed();
		int g=c.getGreen();
		int b=c.getBlue();
		int top=(int) (r*0.3+g*0.59+b*0.11);
		return top;
	}

	private static int setwhite(int rgb){
		Color c=new Color(rgb);
		int r=c.getRed();
		int g=c.getGreen();
		int b=c.getBlue();
		if(r==g&&g==b&&r>224){
			int max;
			return max=new Color(255,255,255).getRGB();
		}
		if(r>126) {
			int max;
			return max=new Color(255,255,255).getRGB();
		}
		return rgb;
	}
}
