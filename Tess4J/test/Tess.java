import java.io.File;  
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.awt.Color;
import java.awt.Graphics2D;
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.imageio.ImageIO;

import LSD.LSD;
import LSD.Line;
import net.sourceforge.tess4j.*;
import net.sourceforge.tess4j.util.ImageHelper;
  
public class Tess {  
	 static String [] retainWord = new String[]{"企业注册号","企业名称","类型","法人","签证机关","校准时间","经营范围","成立时间","注册资本","住所"};//保留字 
     static ArrayList<String> temp1=new   ArrayList<String>();
     static ArrayList<String> temp2=new   ArrayList<String>();
     static boolean tem1=true;
     static boolean tem2=true;
     static boolean tem3=false;
     static boolean tem4=false;
     static  String h=null;
  
     static int depth=1;  
     static int count = 0;
    
  /***遍历指定文件读取读取图片****/    
    public static String find(String pathName,int depth) throws IOException{  
        int filecount=0;  
        //获取pathName的File对象  
        File dirFile = new File(pathName);  
        //判断该文件或目录是否存在，不存在时在控制台输出提醒  
        if (!dirFile.exists()) {  
            System.out.println("do not exit");  
            return "不存在该文件夹";  
        }  
        //判断如果不是一个目录，就判断是不是一个文件，时文件则输出文件路径  
        if (!dirFile.isDirectory()) {  
            if (dirFile.isFile()) {  
                System.out.println(dirFile.getCanonicalFile());  
            }  
            return "不是文件夹";  
        }  
        //获取此目录下的所有文件名与目录名  
        String[] fileList = dirFile.list();  
        int currentDepth=depth+1;  
        for (int i = 0; i < fileList.length; i++) {  
            //遍历文件目录  
            String string = fileList[i];  
            //File("documentName","fileName")是File的另一个构造器  
            File file = new File(dirFile.getPath(),string);  
            String name = file.getName();  
            //如果是一个目录，搜索深度depth++，输出目录名后，进行递归  
            if (file.isDirectory()) {  
                //递归  
                find(file.getCanonicalPath(),currentDepth);  
            }else{  
                //判断是否是图片格式，如果是的话，就提取里面的内容，并将其企业名称和注册号的对应数据存入到相应的列表中
            	if (name.endsWith(".png")||name.endsWith(".jpg")) 
            	{
            		File imageFile = new File(dirFile+"\\"+name);
//                   System.out.println(dirFile+"\\"+name);  
            		BufferedImage img = ImageIO.read(imageFile);
            		//如果是拍摄的照片的话执行该条件--lsz
            		if(name.endsWith(".jpg")) {
            			double max_x1=0;
            			double max_x2=0;
            			double max_y1=0;
            			double max_y2=0;
            			double min_x1=0;
            			double min_x2=0;
            			double min_y1=0;
            			double min_y2=0;

            			double min_x11=0;
            			double min_x22=0;
            			double min_y11=0;
            			double min_y22=0;
            			Graphics2D g2d = img.createGraphics();
            			int x = img.getWidth();
            			int y = img.getHeight();
            			max_x1=x;
            			System.out.println(x);
            			System.out.println(y);
            		
            			HashSet<Line> lines = new HashSet<Line>();

            			
            			double [] arr = img.getData().getPixels(0,0,x,y,new double[x*y*3]);

            			double [] arr2 = new double[x*y];
            		
            			System.out.println("arr.length:"+arr.length);
            			int c=0;
            			for(int i1 = 0; i1 < arr.length-3; i1+=3) {
            				double B = arr[i1];
            				double G = arr[i1+1];
            				double R = arr[i1+2];
            				double level = R * 0.2126 + G * 0.7152 + B * 0.0722;
            				arr2[c++] = level;
            			}

            			LSD lsd = new LSD();

            			double [] out = lsd.lsd(arr2,x,y);

            			for(int i1 = 0; i1 < lsd.n_out; i1++) {
            				for (int j = 0; j < 7; j++)
            				
            				lines.add(new Line(out[7 * i1 + 0], out[7 * i1 + 1],
            						out[7 * i1 + 2], out[7 * i1 + 3]));

            			}
            			
            			
            			int z=0;
            			System.out.println("开始");
            			for (Line l : lines) {
            				if((l.x2-l.x1>x/8)&&(l.y1>min_y11)&&(l.y1<y/2)&&(Math.abs(l.y1-l.y2)<10))
            				{
            					System.out.println("true");
            					min_x11=l.x1;
            					min_x22=l.x2;
            					min_y11=l.y1;
            					min_y22=l.y2;
            				}
            				if(l.x2-l.x1>x/8)
            				{	
            				System.out.println(l.y2-l.y1);
            				System.out.println(l.x1);
            				System.out.println(l.x2);
            				System.out.println(l.y1);
            				System.out.println(l.y2);
            				System.out.println(++z);}
            				if((l.y2-l.y1>y/8)&&(l.x1>min_x1)&&(l.x1<x/2)&&(Math.abs(l.x1-l.x2)<10))
            				{
            					min_x1=l.x1;
            					min_x2=l.x2;
            					min_y1=l.y1;
            					min_y2=l.y2;
            				}
            				else if((l.y2-l.y1>y/8)&&(l.x1<max_x1)&&(l.x1>x/2)&&(Math.abs(l.x1-l.x2)<10))
            				{
            					
            					max_x1=l.x1;
            					max_x2=l.x2;
            					max_y1=l.y1;
            					max_y2=l.y2;
            				}
            				
            			}
            			g2d.setPaint(Color.red);
            		//	g2d.drawLine((int)max_x1,(int)max_y1,(int)max_x2,(int)max_y2);
            		
            		
            			//g2d.drawLine((int)min_x1,(int)min_y1,(int)min_x2,(int)min_y2);
            			g2d.drawLine((int)min_x11,(int)min_y11,(int)min_x22,(int)min_y22);
            		
            			System.out.println(max_x1+" "+min_x1+" "+min_y11);
            			System.out.println(x/2);
            			System.out.println(y/8);
            			System.out.println("结束");
            			System.out.println(lines.size());
            			img=img.getSubimage((int)min_x1,(int)min_y11,(int)(max_x1-min_x1),(int)(y*2/3));
            		}
//            		img = ImageHelper.getScaledInstance(img, img.getHeight() * 2, img.getWidth() * 2);
//            		System.out.println(img);
            		ITesseract instance = new Tesseract();  // JNA Interface Mapping    
//                     ITesseract instance = new Tesseract1(); // JNA Direct Mapping  
            		//代码的文件夹被
            		//原本的请使用你们的路径
            		instance.setDatapath("D:\\EclipseWold\\Tess4J-3.4.6-src\\Tess4J\\tessdata");
            		instance.setLanguage("chi_sim");//添加中文字库   
            		try {    
            			String result =instance.doOCR(img) ;
            			String [] classify=result.split("/n");
            			
            			for(int i1=0;i1<classify.length;i1++)
            			{
            				if(classify[i1].toString().indexOf("企业注册号")!=-1||classify[i1].toString().indexOf("企业名称")!=-1)
            				{
            					if(classify[i1].toString().indexOf("企业名称")!=-1)
            					{   
            						h=classify[i1].substring(classify[i1].indexOf(':')+1 );
            						
            						while(!tem3&&i1<classify.length-1) {
            							i1++;
            						for(int j=0;j<retainWord.length;j++)
            						{	if(classify[i1].toString().indexOf(retainWord[j])!=-1)
            						   {
            							  tem3=true;
            							  break;
            						   }
            						}
            						if(tem3)
            						{
            							tem3=false;
            							i1--;
            						}
            						else
            							h=h+classify[i1];
            						}
            					
            						temp1.add(h);
            						tem1=false;
            					}
            					else
            					{
                                  h=classify[i1].substring(classify[i1].indexOf(':')+1 );
            						temp2.add(classify[i1]);
            						tem2=false;
            					}
            					
            				}
            				
            			}
            			if(tem1)
            				temp1.add("none");
            			if(tem2)
            				temp2.add("none");
            			//	System.out.println(classify[i]); 
            			System.out.println(result);    
            			Main.showMessage.setText(result+"\n"+"已经识别出:"+(++count)+"张");
            		} catch (TesseractException e) {    
            			System.err.println(e.getMessage());    
            		} 
            	}
            }  
        }  
        return "完成";
    }  
      
   /* public static void main(String[] args) throws IOException{  
        find("././b", depth);  //遍历图片，提取指定内容，这里是图片相对位置的存储路径不需要修改
  	    PoiExcel poiexcel=new PoiExcel();
	    poiexcel.poi(temp1,temp2);//生成Excel表格
		StartExcel startexcel=new StartExcel();
	    startexcel.start();//自动打开Excel表格
    }  */
}  