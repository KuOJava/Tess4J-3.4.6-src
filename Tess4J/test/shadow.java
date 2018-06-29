import java.awt.Color;
import java.awt.image.BufferedImage;
import java.io.IOException;
import java.util.Vector;
import net.sourceforge.tess4j.ITesseract;
import net.sourceforge.tess4j.Tesseract;
import net.sourceforge.tess4j.TesseractException;

//import testocr.shadow.peek_range;

public class shadow {
	
	static class peek_range{
		int begin;
		int end;
	}
	
	public String ocr(BufferedImage img) throws IOException {
		Invert inv=new Invert();
		img=Invert.invert(img);
		HoughLine HL = new HoughLine();
        img=HL.hough(img);
		String results="";
		Vector<Integer> a=new Vector<>();
		Vector<peek_range> peek=new Vector<>();
		horizational(img,a);
		//这里
		GetPeekRange(peek,a,1,8);
		BufferedImage[] nbi = new BufferedImage[peek.size()];
		Vector<Integer>[] b=new Vector[peek.size()];
		for(int i=0;i<peek.size();i++){
			int end=peek.get(i).end,begin=peek.get(i).begin;
			nbi[i]=new BufferedImage(img.getWidth(),end-begin+1,BufferedImage.TYPE_BYTE_BINARY);
			for(int j=0;j<img.getWidth();j++)
				for(int k=begin;k<=end;k++){
					nbi[i].setRGB(j, k-begin, img.getRGB(j, k));
				}
			b[i]=new Vector<Integer>();
			vertical(nbi[i],b[i]);
		}
		//这里
		int num[]=new int[30];
		BufferedImage im = null;
		Vector< Vector<peek_range> > peek_v=new Vector<>();
		Vector< peek_range > peek_v_x=new Vector<>();
		for(int i=0;i<peek.size();i++){
			GetPeekRange(peek_v_x,b[i],1,1);
			peek_v.insertElementAt((Vector< peek_range >)peek_v_x.clone(), i);
			peek_v_x.clear();
			for(int j=0;j<peek_v.elementAt(i).size();j++){
				int end=peek_v.elementAt(i).get(j).end,begin=peek_v.elementAt(i).get(j).begin;
				if(end-begin+1<30)
					num[end-begin+1]++;
			}
		}
		
		int max=0,average=0;
		for(int i=0;i<30;i++)
			if(num[i]>max)
				average=i;
		for(int i=0;i<peek.size()*0.4;i++){
			for(int j=0;j<peek_v.elementAt(i).size()-1;j++){
				int end=peek_v.elementAt(i).get(j+1).end,begin=peek_v.elementAt(i).get(j).begin;
				if(end-begin+1<average){
					peek_range p=new peek_range();
					p.begin = begin;
					p.end = end;
					peek_v.elementAt(i).set(j, p);
					peek_v.elementAt(i).remove(j+1);
				}
			}
		}
		ITesseract instance = new Tesseract();  // JNA Interface Mapping    
		instance.setDatapath("./Tess4J/tessdata");
		instance.setLanguage("chi12");//添加中文字库   
		
		for(int i=0;i<peek.size()*0.5;i++){
			for(int j=0;j<peek_v.elementAt(i).size();j++){
				int end=peek_v.elementAt(i).get(j).end,begin=peek_v.elementAt(i).get(j).begin;
				if((end-begin+1)*1.0/nbi[i].getHeight()>1.5){
					int newEnd=begin+average-2;
					peek_range p=new peek_range(),q=new peek_range();
					p.begin = begin;
					q.begin=p.end = newEnd;
					peek_v.elementAt(i).set(j, p);
					q.end=end;
					peek_v.elementAt(i).add(j+1,q);
				}
			}
			if(peek_v.elementAt(i).size()>0){
				im=new BufferedImage(peek_v.elementAt(i).get(peek_v.elementAt(i).size()-1).end-peek_v.elementAt(i).get(0).begin
						+3*peek_v.elementAt(i).size(),nbi[i].getHeight(),BufferedImage.TYPE_BYTE_BINARY);
				for(int j=0;j<peek_v.elementAt(i).size();j++){
					int end=peek_v.elementAt(i).get(j).end,begin=peek_v.elementAt(i).get(j).begin;
					for(int k=begin;k<=end;k++){
						for(int s=0;s<nbi[i].getHeight();s++)
							im.setRGB(k+j*2-peek_v.elementAt(i).get(0).begin,s, nbi[i].getRGB(k,s));
					}				
				}
				im=inv.invert(im);
				try {    
					results+=instance.doOCR(im);					
				} catch (TesseractException e) {    
					System.err.println(e.getMessage());    
				}
			}						
		}
		
		return results;
		
	}	

	public static void horizational(BufferedImage img,Vector<Integer> a){
		int num=0,max=new Color(255,255,255).getRGB();
		for(int i=0;i<img.getHeight();i++){
			for(int j=0;j<img.getWidth();j++){
				if(img.getRGB(j,i)==max)
					num++;
			}
			a.add(num);
			num=0;
		}
	}
	
	public static void vertical(BufferedImage img,Vector<Integer> b){
		int num=0,max=new Color(255,255,255).getRGB();
		for(int i=0;i<img.getWidth();i++){
			for(int j=0;j<img.getHeight();j++){
				if(img.getRGB(i,j)==max)
					num++;
			}
			b.add(num);
			num=0;
		}
	}
	
	public static void GetPeekRange(Vector<peek_range> peek,Vector a,int min_thresh,int  min_range ){
		  int begin = 0;
		  int end = 0;
		  for (int i = 0; i < a.size(); i++){
			  if ((int)a.get(i) > min_thresh && begin == 0){
		     	begin = i;
			  }
			  else if ((int)a.get(i)> min_thresh && begin != 0){
		            continue;
		      }
			  else if ((int)a.get(i) < min_thresh && begin != 0){
				  end = i;
				  if (end - begin >= min_range){
					  peek_range p=new peek_range();
					  p.begin = begin;
					  p.end = end;
					  peek.add(p);
					  begin = 0;
				  }
			  }
			  else if ((int)a.get(i) < min_thresh || begin == 0){
		            continue;
			  }
		  }
	}
}