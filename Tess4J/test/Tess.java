import java.io.File;  
import java.io.IOException;
import java.util.ArrayList;  
import java.awt.image.BufferedImage;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import javax.imageio.ImageIO;
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
  
    private static int depth=1;  
    
    
  /***遍历指定文件读取读取图片****/    
    public static void find(String pathName,int depth) throws IOException{  
        int filecount=0;  
        //获取pathName的File对象  
        File dirFile = new File(pathName);  
        //判断该文件或目录是否存在，不存在时在控制台输出提醒  
        if (!dirFile.exists()) {  
            System.out.println("do not exit");  
            return ;  
        }  
        //判断如果不是一个目录，就判断是不是一个文件，时文件则输出文件路径  
        if (!dirFile.isDirectory()) {  
            if (dirFile.isFile()) {  
                System.out.println(dirFile.getCanonicalFile());  
            }  
            return ;  
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
//            		img = ImageHelper.getScaledInstance(img, img.getHeight() * 2, img.getWidth() * 2);
//            		System.out.println(img);
            		ITesseract instance = new Tesseract();  // JNA Interface Mapping    
//                     ITesseract instance = new Tesseract1(); // JNA Direct Mapping  
            		//代码的文件夹被修改
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
            		} catch (TesseractException e) {    
            			System.err.println(e.getMessage());    
            		} 
            	}
            }  
        }  
    }  
      
    public static void main(String[] args) throws IOException{  
        find("././b", depth);  //遍历图片，提取指定内容，这里是图片相对位置的存储路径不需要修改
  	    PoiExcel poiexcel=new PoiExcel();
	    poiexcel.poi(temp1,temp2);//生成Excel表格
		StartExcel startexcel=new StartExcel();
	    startexcel.start();//自动打开Excel表格
    }  
}  