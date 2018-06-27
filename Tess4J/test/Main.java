import java.awt.Component;
import java.awt.Container;
import java.awt.GridBagConstraints;
import java.awt.GridBagLayout;
import java.awt.event.ActionEvent;
import java.awt.event.ActionListener;
import java.io.IOException;

import javax.swing.JButton;
import javax.swing.JFileChooser;
import javax.swing.JFrame;
import javax.swing.JLabel;
import javax.swing.JTextArea;
import javax.swing.JTextField;

//这是主函数，图形用户界面编写在这里
public class Main {
	JTextField inputPath;
	private JButton choosePath;
	private JButton confirm;
	private JLabel situation;
	public static JTextArea showMessage;
	private JLabel jumpToExcel;
	private JButton clickToCopy;

	public Main() {
		// 创建一个框架作为顶层容器
		// JFrame.setDefaultLookAndFeelDecorated(true);
		JFrame f = new JFrame("工商图片识别系统");
		// f.setDefaultCloseOperation(JFrame.EXIT_ON_CLOSE);
		// 获取框架中的内容面板
		Container contentPane = f.getContentPane();
		// 设置网格布局
		contentPane.setLayout(new GridBagLayout());
		// 使用网格布局约束器
		GridBagConstraints c = new GridBagConstraints();
		// 界面设计,根据需求文档中提供的图形界面得
		// '输入路径'文本框：处于第1行，占3列;'选择路径'按钮：处于第1行，第4列；'确定'按钮：处于第1行，第5列
		// '识别情况'提示框：处于第2行，占所有列
		// 显示消息的文本框：处于第3,4行，占所有列
		// '点击进入输出的Excel的根目录'超链接：处于第5行，占3列；'点击复制文件'按钮：处于第5行，占两列

		// '输入路径'的文本框
		inputPath = new JTextField(1);
		// 第一行第一列开始
		c.gridx = 0;
		c.gridy = 0;
		// 跨3列
		c.gridwidth = 3;
		c.gridheight = 1;
		c.fill = GridBagConstraints.BOTH;
		c.weightx = 3;
		c.weighty = 1;
		// 设置点击事件

		// 根据c的设置添加输入框
		contentPane.add(inputPath, c);

		// '选择路径'按钮
		choosePath = new JButton("选择路径");
		// 选择路径点击事件
		choosePath.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				JFileChooser chooser = new JFileChooser();
				chooser.setFileSelectionMode(JFileChooser.DIRECTORIES_ONLY);
				Component parent = null;
				int returnVal = chooser.showOpenDialog(parent);
				if (returnVal == JFileChooser.APPROVE_OPTION) {
					inputPath.setText(chooser.getSelectedFile().getAbsolutePath());
				}
			}
		});
		// 放在第一行第四列
		c.gridx = 3;
		c.gridy = 0;
		c.gridwidth = 1;
		c.weightx = 1;
		c.weighty = 1;
		// c.gridheight=1;----不用赋值
		contentPane.add(choosePath, c);
		// '确定'按钮
		confirm = new JButton("确定");
		// 确定按钮的点击事件
		confirm.addActionListener(new ActionListener() {

			@Override
			public void actionPerformed(ActionEvent e) {
				// 点击后获取路径读取内容
					String path = inputPath.getText();
					if (path.equals("")) {
						showMessage.setText("你没有输入任何路径\n");
					} else {
						showMessage.setText("正在识别。。。。。。");
						new Thread(new Runnable() {

							@Override
							public void run() {
								// TODO Auto-generated method stub
								String message="";
								try {
									message = Tess.find(path, Tess.depth);
									
								} catch (IOException e) {
									// TODO Auto-generated catch block
									e.printStackTrace();
								}
								PoiExcel poiexcel = new PoiExcel();
								poiexcel.poi(Tess.temp1, Tess.temp2);// 生成Excel表格
								showMessage.setText(message);
							}
						}).start();

					}
			}
		});
		c.gridx = 4;
		// c.gridy=4;----不用赋值
		// c.gridwidth=1;----不用赋值
		// c.gridheight=1;----不用赋值
		c.weightx = 1;
		c.weighty = 1;
		contentPane.add(confirm, c);

		// '识别情况:'提示框
		situation = new JLabel("识别情况:");
		c.gridx = 0;
		c.gridy = 1;
		c.weightx = 1;
		c.weighty = 1;
		// c.gridwidth=1;----不用赋值
		// c.gridheight=1;----不用赋值
		contentPane.add(situation, c);

		// '显示消息的框'
		showMessage = new JTextArea("");
		showMessage.setSize(200, 200);
		// c.gridx=0;
		c.gridy = 2;
		c.gridwidth = 5;
		c.gridheight = 2;
		c.weightx = 5;
		c.weighty = 2;
		contentPane.add(showMessage, c);

		// '点击进入输出的Excel的根目录'
		jumpToExcel = new JLabel("点击进入输出的Excel的根目录");
		c.gridx = 0;
		c.gridy = 4;
		c.gridwidth = 3;
		c.gridheight = 1;
		c.weightx = 3;
		c.weighty = 1;
		contentPane.add(jumpToExcel, c);

		// '点击复制文件'
		clickToCopy = new JButton("点击复制文件");
		c.gridx = 3;
		c.gridwidth = 2;
		c.weightx = 2;
		c.weighty = 1;
		contentPane.add(clickToCopy, c);
		// 显示窗口
		f.pack();
		f.setSize(800, 1000);
		f.setVisible(true);
	}

	public static void main(String[] args) {
		new Main();
	}
}
