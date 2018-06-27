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
	 static String [] retainWord = new String[]{"��ҵע���","��ҵ����","����","����","ǩ֤����","У׼ʱ��","��Ӫ��Χ","����ʱ��","ע���ʱ�","ס��"};//������ 
     static ArrayList<String> temp1=new   ArrayList<String>();
     static ArrayList<String> temp2=new   ArrayList<String>();
     static boolean tem1=true;
     static boolean tem2=true;
     static boolean tem3=false;
     static boolean tem4=false;
     static  String h=null;
  
     static int depth=1;  
    
    
  /***����ָ���ļ���ȡ��ȡͼƬ****/    
    public static void find(String pathName,int depth) throws IOException{  
        int filecount=0;  
        //��ȡpathName��File����  
        File dirFile = new File(pathName);  
        //�жϸ��ļ���Ŀ¼�Ƿ���ڣ�������ʱ�ڿ���̨�������  
        if (!dirFile.exists()) {  
            System.out.println("do not exit");  
            return ;  
        }  
        //�ж��������һ��Ŀ¼�����ж��ǲ���һ���ļ���ʱ�ļ�������ļ�·��  
        if (!dirFile.isDirectory()) {  
            if (dirFile.isFile()) {  
                System.out.println(dirFile.getCanonicalFile());  
            }  
            return ;  
        }  
        //��ȡ��Ŀ¼�µ������ļ�����Ŀ¼��  
        String[] fileList = dirFile.list();  
        int currentDepth=depth+1;  
        for (int i = 0; i < fileList.length; i++) {  
            //�����ļ�Ŀ¼  
            String string = fileList[i];  
            //File("documentName","fileName")��File����һ��������  
            File file = new File(dirFile.getPath(),string);  
            String name = file.getName();  
            //�����һ��Ŀ¼���������depth++�����Ŀ¼���󣬽��еݹ�  
            if (file.isDirectory()) {  
                //�ݹ�  
                find(file.getCanonicalPath(),currentDepth);  
            }else{  
                //�ж��Ƿ���ͼƬ��ʽ������ǵĻ�������ȡ��������ݣ���������ҵ���ƺ�ע��ŵĶ�Ӧ���ݴ��뵽��Ӧ���б���
            	if (name.endsWith(".png")||name.endsWith(".jpg")) 
            	{
            		File imageFile = new File(dirFile+"\\"+name);
//                   System.out.println(dirFile+"\\"+name);  
            		BufferedImage img = ImageIO.read(imageFile);
//            		img = ImageHelper.getScaledInstance(img, img.getHeight() * 2, img.getWidth() * 2);
//            		System.out.println(img);
            		ITesseract instance = new Tesseract();  // JNA Interface Mapping    
//                     ITesseract instance = new Tesseract1(); // JNA Direct Mapping  
            		//������ļ��б�
            		//ԭ������ʹ�����ǵ�·��
            		instance.setDatapath("D:\\EclipseWold\\Tess4J-3.4.6-src\\Tess4J\\tessdata");
            		instance.setLanguage("chi_sim");//��������ֿ�   
            		try {    
            			String result =instance.doOCR(img) ;
            			String [] classify=result.split("/n");
            			
            			for(int i1=0;i1<classify.length;i1++)
            			{
            				if(classify[i1].toString().indexOf("��ҵע���")!=-1||classify[i1].toString().indexOf("��ҵ����")!=-1)
            				{
            					if(classify[i1].toString().indexOf("��ҵ����")!=-1)
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
      
   /* public static void main(String[] args) throws IOException{  
        find("././b", depth);  //����ͼƬ����ȡָ�����ݣ�������ͼƬ���λ�õĴ洢·������Ҫ�޸�
  	    PoiExcel poiexcel=new PoiExcel();
	    poiexcel.poi(temp1,temp2);//����Excel���
		StartExcel startexcel=new StartExcel();
	    startexcel.start();//�Զ���Excel���
    }  */
}  